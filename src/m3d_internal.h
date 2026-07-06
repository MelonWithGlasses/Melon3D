// Melon3D internal data structures. Not part of the public API.
//
// Architecture summary:
//  - broadphase: sorted spatial hash grid, rebuilt every step (no BVH)
//  - narrowphase: fresh contact manifolds every step, speculative margin
//    scaled by relative velocity (no persistent warm-started manifolds)
//  - solver: XPBD (extended position based dynamics) with substeps -
//    positional projection of contacts/joints, velocity reconstruction
//    from positions, then a velocity pass for restitution and dynamic
//    friction (Mueller et al., "Detailed Rigid Body Simulation with XPBD")
#pragma once

#include "melon3d.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace m3d
{

// Tuning constants
constexpr float kLinearSlop = 0.005f;
constexpr float kSpeculativeBase = 4.0f * kLinearSlop;
constexpr float kMaxDepenetrationSpeed = 3.0f; // m/s, cap on positional projection
// Static friction position corrections are under-relaxed. A full correction
// acts like a base re-centering controller with loop gain > 1: tall towers
// pump their sway mode and never sleep. The velocity-pass dynamic friction
// covers the remainder.
constexpr float kFrictionRelaxation = 0.5f;
constexpr float kAABBMargin = 0.1f;
constexpr float kTimeToSleep = 0.5f;
constexpr float kSleepLinearVelSq = 0.05f * 0.05f;  // (m/s)^2
constexpr float kSleepAngularVelSq = 0.05f * 0.05f; // (rad/s)^2
// Sleep decisions use a low-pass-filtered energy |v|^2 + |w|^2 with
// hysteresis instead of instantaneous velocities: XPBD velocities are
// reconstructed from positions and carry single-step noise that would
// otherwise keep resetting the sleep timers of settling piles.
constexpr float kSleepEnergyThresh = kSleepLinearVelSq + kSleepAngularVelSq;
constexpr float kSleepEnergyAlpha = 0.25f; // per-step smoothing factor
constexpr float kSleepWakeFactor = 2.0f;   // reset timers above this multiple
constexpr int32_t kOversizedCellSpan = 64;          // grid cells before a body goes to the coarse list
// Islands with at least this many contacts are solved with graph coloring:
// contacts are bucketed so no two in a bucket share a dynamic body, giving
// parallelism INSIDE one island (a pyramid is a single island). Few colors
// keep buckets large and the per-color barrier count low; the leftovers go
// to a serial overflow pass.
constexpr int32_t kBigIslandContacts = 300;
constexpr int32_t kMaxGraphColors = 8;

// ---------------------------------------------------------------------------
// Bodies
// ---------------------------------------------------------------------------

struct Shape
{
	m3d_shape_type type;
	float radius;        // sphere, capsule
	float halfHeight;    // capsule (local Y)
	m3d_vec3 halfExtents; // box
};

struct Body
{
	m3d_vec3 position; // center of mass (shapes are centered on the body origin)
	m3d_quat rotation;
	m3d_vec3 linearVelocity;
	m3d_vec3 angularVelocity;

	// XPBD substep state
	m3d_vec3 prevPosition;
	m3d_quat prevRotation;
	m3d_vec3 prevLinearVelocity;
	m3d_vec3 prevAngularVelocity;

	float invMass;
	m3d_vec3 invInertiaLocal;  // diagonal in body frame
	m3d_mat3 invInertiaWorld;  // refreshed during solve

	Shape shape;
	float friction;
	float restitution;
	float mass;

	m3d_body_type type;
	float linearDamping;
	float angularDamping;
	float gravityScale;

	float sleepTime;
	float smoothedEnergy; // low-pass |v|^2 + |w|^2 for sleep decisions
	bool enableSleep;
	bool isAwake;
	bool inUse;
	uint32_t generation;
	uint64_t userData;
};

inline m3d_transform BodyTransform(const Body& b)
{
	return { b.position, b.rotation };
}

inline void UpdateInvInertiaWorld(Body& b)
{
	// I^-1_world = R * diag(invInertiaLocal) * R^T
	m3d_mat3 r = m3d_quat_to_mat3(b.rotation);
	m3d_vec3 d = b.invInertiaLocal;
	m3d_mat3 rd;
	rd.cx = m3d_scale(d.x, r.cx);
	rd.cy = m3d_scale(d.y, r.cy);
	rd.cz = m3d_scale(d.z, r.cz);
	m3d_mat3 m;
	m.cx = m3d_v3(rd.cx.x * r.cx.x + rd.cy.x * r.cy.x + rd.cz.x * r.cz.x,
				  rd.cx.y * r.cx.x + rd.cy.y * r.cy.x + rd.cz.y * r.cz.x,
				  rd.cx.z * r.cx.x + rd.cy.z * r.cy.x + rd.cz.z * r.cz.x);
	m.cy = m3d_v3(rd.cx.x * r.cx.y + rd.cy.x * r.cy.y + rd.cz.x * r.cz.y,
				  rd.cx.y * r.cx.y + rd.cy.y * r.cy.y + rd.cz.y * r.cz.y,
				  rd.cx.z * r.cx.y + rd.cy.z * r.cy.y + rd.cz.z * r.cz.y);
	m.cz = m3d_v3(rd.cx.x * r.cx.z + rd.cy.x * r.cy.z + rd.cz.x * r.cz.z,
				  rd.cx.y * r.cx.z + rd.cy.y * r.cy.z + rd.cz.y * r.cz.z,
				  rd.cx.z * r.cx.z + rd.cy.z * r.cy.z + rd.cz.z * r.cz.z);
	b.invInertiaWorld = m;
}

m3d_aabb ComputeShapeAABB(const Shape& shape, m3d_transform xf);

// ---------------------------------------------------------------------------
// Contacts (regenerated every step; only touch state persists for events)
// ---------------------------------------------------------------------------

struct ManifoldPoint
{
	// anchors in each body's local frame (relative to center of mass)
	m3d_vec3 localAnchorA;
	m3d_vec3 localAnchorB;
	// stiction targets for static friction. Inherited across manifold
	// regenerations (matched by feature id) while the pair has not slipped,
	// so the friction target never ratchets forward with micro-creep.
	m3d_vec3 frictionAnchorA;
	m3d_vec3 frictionAnchorB;
	float separation; // at detection; < 0 means penetration
	uint32_t id;

	// per-substep XPBD state
	float lambdaN; // accumulated positional normal correction
	float lambdaT; // accumulated positional friction correction
	float separationOffset;

	// frozen at the start of each substep so the projection iterations and
	// the velocity pass avoid re-rotating anchors and recomputing masses
	m3d_vec3 rA; // world anchor offset from body A center
	m3d_vec3 rB;
	float wN; // generalized inverse mass along the normal
	// frozen angular response Jacobians J = I^-1 (r x n): applying a normal
	// correction dLambda changes the rotation by 0.5*(dLambda*J) x q without
	// touching the inertia matrix or cross products in the hot loop
	m3d_vec3 angJA;
	m3d_vec3 angJB;
};

struct Manifold
{
	m3d_vec3 normal; // world, points from A to B
	int32_t pointCount;
	ManifoldPoint points[4];
};

struct Contact
{
	int32_t bodyA; // bodyA < bodyB
	int32_t bodyB;
	uint64_t key;
	Manifold manifold;
	float friction;
	float restitution;
	bool touching;

	// Pose snapshot at manifold generation. While neither body has moved
	// beyond a small tolerance the manifold is reused: this both skips
	// narrowphase for resting contacts and turns the anchors into persistent
	// stiction targets, so static friction converges to zero drift instead
	// of accumulating creep.
	bool hasCache;
	m3d_vec3 cachePosA, cachePosB;
	m3d_quat cacheRotA, cacheRotB;
};

inline bool PoseChanged(const Body& b, m3d_vec3 cachePos, m3d_quat cacheRot)
{
	constexpr float kTolPosSq = 1e-3f * 1e-3f;      // 1 mm
	constexpr float kTolHalfAngleSq = 0.25f * 1e-6f; // ~1 mrad
	if (m3d_length_sq(m3d_sub(b.position, cachePos)) > kTolPosSq)
	{
		return true;
	}
	m3d_quat d = m3d_quat_mul(b.rotation, m3d_quat_conj(cacheRot));
	return d.x * d.x + d.y * d.y + d.z * d.z > kTolHalfAngleSq;
}

inline uint64_t ContactKey(int32_t a, int32_t b)
{
	return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

// index of the lowest set bit; caller guarantees v != 0
inline int32_t CountTrailingZeros32(uint32_t v)
{
#if defined(_MSC_VER)
	unsigned long idx;
	_BitScanForward(&idx, v);
	return (int32_t)idx;
#else
	return __builtin_ctz(v);
#endif
}

// Compute the contact manifold between two shapes. Normal points from A to B.
// Points with separation up to `margin` are kept as speculative contacts.
void CollideShapes(const Shape& shapeA, m3d_transform xfA, const Shape& shapeB, m3d_transform xfB, float margin,
				   Manifold* manifold);

// Carry stiction friction anchors from the previous manifold into the fresh
// one (matched by feature id) unless the pair slipped too far.
void InheritFrictionAnchors(const Manifold& oldManifold, Manifold& manifold, const Body& a, const Body& b);

// Ray cast against a single shape. translation = full ray extent, fraction in [0,1].
bool RayCastShape(const Shape& shape, m3d_transform xf, m3d_vec3 origin, m3d_vec3 translation, float* fraction,
				  m3d_vec3* normal);

// ---------------------------------------------------------------------------
// Joints (solved positionally, XPBD compliance)
// ---------------------------------------------------------------------------

enum class JointType
{
	Distance,
	Ball,
};

struct Joint
{
	JointType type;
	int32_t bodyA;
	int32_t bodyB;
	m3d_vec3 localAnchorA;
	m3d_vec3 localAnchorB;

	// distance
	float length;
	float hertz;
	float dampingRatio;

	// transient
	float compliance; // XPBD alpha [m/N], 0 = rigid
	float lambda;     // per-substep accumulator

	bool inUse;
	uint32_t generation;
};

// ---------------------------------------------------------------------------
// Broadphase: sorted spatial hash grid
// ---------------------------------------------------------------------------

// Two-tier sorted spatial hash. Awake bodies go to the ACTIVE tier, rebuilt
// every step. Static and sleeping bodies go to the STABLE tier, rebuilt only
// when membership changes. Pair generation never enumerates stable x stable:
// sleeping piles cost nothing in the broadphase.
class SpatialGrid
{
public:
	struct Entry
	{
		uint64_t key;
		int32_t body;
	};

	void BuildActive(float cellSize, const std::vector<m3d_aabb>& aabbs, const std::vector<uint8_t>& isActive);
	void BuildStable(float cellSize, const std::vector<m3d_aabb>& aabbs, const std::vector<uint8_t>& isStable);

	// callback(bodyA, bodyB) for every candidate pair with at least one
	// active body (deduplicated)
	template <typename F>
	void QueryPairs(const std::vector<m3d_aabb>& aabbs, F&& callback) const
	{
		size_t n = m_active.size();
		size_t runStart = 0;
		while (runStart < n)
		{
			size_t runEnd = runStart + 1;
			uint64_t key = m_active[runStart].key;
			while (runEnd < n && m_active[runEnd].key == key)
			{
				++runEnd;
			}

			// active x active within this cell
			for (size_t i = runStart; i < runEnd; ++i)
			{
				for (size_t j = i + 1; j < runEnd; ++j)
				{
					int32_t a = m_active[i].body;
					int32_t b = m_active[j].body;
					// dedupe: a pair sharing k cells is only reported from the
					// cell containing the min corner of the AABB intersection
					m3d_vec3 lo = m3d_max(aabbs[a].lowerBound, aabbs[b].lowerBound);
					if (CellKeyOf(lo) == key)
					{
						callback(a, b);
					}
				}
			}

			// active x stable for the same cell (binary search the stable tier)
			auto it = std::lower_bound(m_stable.begin(), m_stable.end(), key,
									   [](const Entry& e, uint64_t k) { return e.key < k; });
			for (; it != m_stable.end() && it->key == key; ++it)
			{
				for (size_t i = runStart; i < runEnd; ++i)
				{
					int32_t a = m_active[i].body;
					int32_t b = it->body;
					if (a == b)
					{
						continue; // stale stable entry of a body that woke
					}
					m3d_vec3 lo = m3d_max(aabbs[a].lowerBound, aabbs[b].lowerBound);
					if (CellKeyOf(lo) == key)
					{
						callback(a, b);
					}
				}
			}
			runStart = runEnd;
		}

		// oversized bodies are tested against everything relevant
		for (size_t oi = 0; oi < m_activeOversized.size(); ++oi)
		{
			int32_t a = m_activeOversized[oi];
			for (int32_t b = 0; b < (int32_t)aabbs.size(); ++b)
			{
				if (b != a && (m_inActive[b] || m_inStable[b]))
				{
					callback(a, b);
				}
			}
			for (size_t oj = oi + 1; oj < m_activeOversized.size(); ++oj)
			{
				callback(a, m_activeOversized[oj]);
			}
			for (int32_t b : m_stableOversized)
			{
				if (b != a)
				{
					callback(a, b);
				}
			}
		}
		for (int32_t a : m_stableOversized)
		{
			for (size_t i = 0; i < m_inActive.size(); ++i)
			{
				if (m_inActive[i])
				{
					callback(a, (int32_t)i);
				}
			}
		}
	}

	uint64_t CellKeyOf(m3d_vec3 p) const;
	const std::vector<Entry>& ActiveEntries() const { return m_active; }
	const std::vector<Entry>& StableEntries() const { return m_stable; }
	const std::vector<int32_t>& ActiveOversized() const { return m_activeOversized; }
	const std::vector<int32_t>& StableOversized() const { return m_stableOversized; }
	float CellSize() const { return m_cellSize; }

private:
	static void BuildTier(float cellSize, const std::vector<m3d_aabb>& aabbs, const std::vector<uint8_t>& include,
						  std::vector<Entry>& entries, std::vector<int32_t>& oversized, std::vector<bool>& inGrid);

	float m_cellSize = 1.5f;
	std::vector<Entry> m_active; // sorted by key
	std::vector<Entry> m_stable; // sorted by key
	std::vector<int32_t> m_activeOversized;
	std::vector<int32_t> m_stableOversized;
	std::vector<bool> m_inActive; // in active cells (not oversized)
	std::vector<bool> m_inStable;
};

// ---------------------------------------------------------------------------
// SIMD contact packets (colored solve). 8 contacts = 8 lanes; manifold
// points are iterated sequentially inside the packet (points of one contact
// share its bodies), while the 8 lanes are body-disjoint by coloring. All
// lane math is written as plain float[8] loops that GCC vectorizes fully
// with -mavx2 (branchless: masks are 0/1 multipliers).
// ---------------------------------------------------------------------------

struct alignas(32) ContactPacket8
{
	int32_t contactIndex[8]; // -1 = padding lane
	int32_t bodyA[8];
	int32_t bodyB[8];
	int32_t lanePoints[8];
	int32_t maxPoints;

	float invMassA[8];
	float invMassB[8];
	float nx[8], ny[8], nz[8];

	// per manifold point j (up to 4)
	float active[4][8]; // 1 = real point, 0 = padding
	float lax[4][8], lay[4][8], laz[4][8]; // local anchor A
	float lbx[4][8], lby[4][8], lbz[4][8]; // local anchor B
	float sepOff[4][8];
	float invW[4][8]; // 1 / wN (1 for padding)
	float jAx[4][8], jAy[4][8], jAz[4][8]; // angular Jacobian A
	float jBx[4][8], jBy[4][8], jBz[4][8]; // angular Jacobian B
	float lamN[4][8];
};

// ---------------------------------------------------------------------------
// Thread pool
// ---------------------------------------------------------------------------

// Fork-join pool with lock-free job publication and spin-before-sleep
// workers. The graph-colored solver issues ~140 small parallel stages per
// step; a condition-variable handoff (~10us) would eat the gains, so jobs
// are published with a release store of the job id and hot workers pick
// them up from a bounded spin loop (~1us per stage barrier).
class ThreadPool
{
public:
	explicit ThreadPool(int32_t workerCount);
	~ThreadPool();

	void ParallelFor(int32_t count, int32_t grain, const std::function<void(int32_t)>& fn);
	int32_t WorkerCount() const { return (int32_t)m_threads.size() + 1; }

private:
	void WorkerLoop();
	void RunChunks();

	std::vector<std::thread> m_threads;
	std::mutex m_mutex;
	std::condition_variable m_cvWork;

	// job slot: stable while m_jobOpen; workers enter under the mutex, so a
	// straggler can never observe a half-published or retired job
	const std::function<void(int32_t)>* m_fn = nullptr;
	int32_t m_count = 0;
	int32_t m_grain = 1;
	bool m_jobOpen = false; // guarded by m_mutex
	std::atomic<int32_t> m_next{ 0 };

	std::atomic<uint64_t> m_jobId{ 0 }; // bumped on publish; spun on for discovery
	std::atomic<int32_t> m_activeWorkers{ 0 };
	std::atomic<int32_t> m_sleepers{ 0 };
	std::atomic<bool> m_shutdown{ false };
};

} // namespace m3d

// ---------------------------------------------------------------------------
// World
// ---------------------------------------------------------------------------

struct m3d_world
{
	m3d_world_def def;

	std::vector<m3d::Body> bodies;
	std::vector<int32_t> freeBodies;

	std::vector<m3d::Joint> joints;
	std::vector<int32_t> freeJoints;

	std::vector<m3d::Contact> contacts;
	std::unordered_map<uint64_t, int32_t> contactMap;

	m3d::SpatialGrid grid;
	bool gridDirty = true;
	bool stableDirty = true; // stable tier membership changed (sleep/wake/create/destroy)
	bool stableForce = false; // must rebuild now (stale entries would dangle)
	int32_t stepsSinceStableBuild = 1000;
	std::vector<m3d_aabb> expandedAABBs; // persistent fat AABBs (margin + velocity), indexed by body
	// uint8 (not vector<bool>): written concurrently from the parallel AABB loop
	std::vector<uint8_t> isActiveBody;
	std::vector<uint8_t> isStableBody;
	std::vector<uint8_t> hasFatAABB; // fat AABB valid for this body slot
	std::vector<uint8_t> movedFlag;  // breached its fat AABB this step -> re-pair
	int32_t lastNpRegenCount = 1 << 30; // manifold regens last step; gates parallel narrowphase
	std::vector<uint8_t> inStableTier; // migrated into the stable tier at its last build
	std::vector<uint32_t> rayStamps; // per-body visited stamp for ray casts
	uint32_t rayStampCounter = 0;

	std::vector<m3d_contact_event> beginEvents;
	std::vector<m3d_contact_event> endEvents;

	m3d_profile profile;

	m3d::ThreadPool* pool;

	// solver scratch
	std::vector<std::vector<int32_t>> islandContacts;
	std::vector<std::vector<int32_t>> islandJoints;
	std::vector<std::vector<int32_t>> islandBodies;

	// graph coloring scratch (big islands are colored one at a time)
	std::vector<uint32_t> colorMask; // per body: colors already used
	std::vector<std::vector<int32_t>> colorBuckets;
	std::vector<int32_t> overflowContacts;
	std::vector<std::vector<m3d::ContactPacket8>> colorPackets;
};

namespace m3d
{
void SolveIslandXPBD(m3d_world* world, const std::vector<int32_t>& bodyIndices,
					 const std::vector<int32_t>& contactIndices, const std::vector<int32_t>& jointIndices, float dt,
					 int substepCount);

// Graph-colored variant for big islands: stages run as parallel-for over
// color buckets. Within a bucket bodies are disjoint, so results are
// bit-identical for any thread count (including inline single-thread).
void SolveIslandColoredXPBD(m3d_world* world, const std::vector<int32_t>& bodyIndices,
							const std::vector<int32_t>& contactIndices, const std::vector<int32_t>& jointIndices,
							float dt, int substepCount);
}
