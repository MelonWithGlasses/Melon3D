// World management and the step pipeline:
//   rebuild spatial hash grid -> candidate pairs -> parallel narrowphase ->
//   island building (union-find) -> parallel XPBD island solve -> sleeping.
#include "m3d_internal.h"

#include <algorithm>
#include <chrono>
#include <float.h>
#include <string.h>
#include <unordered_set>

using namespace m3d;

namespace
{

using Clock = std::chrono::steady_clock;

float ToMs(Clock::time_point t0, Clock::time_point t1)
{
	return std::chrono::duration<float, std::milli>(t1 - t0).count();
}

Body* GetBodyChecked(m3d_world* world, m3d_body_id id)
{
	if (id.index < 0 || id.index >= (int32_t)world->bodies.size())
	{
		return nullptr;
	}
	Body& b = world->bodies[id.index];
	if (!b.inUse || b.generation != id.generation)
	{
		return nullptr;
	}
	return &b;
}

m3d_body_id MakeBodyId(m3d_world* world, int32_t index)
{
	return { index, world->bodies[index].generation };
}

void ComputeMassProperties(Body& body, const m3d_shape_def* def)
{
	float density = def->density > 0.0f ? def->density : 1000.0f;
	float mass = 0.0f;
	m3d_vec3 inertia = m3d_v3_zero();

	switch (def->type)
	{
		case M3D_SHAPE_SPHERE:
		{
			float r = def->radius;
			mass = density * (4.0f / 3.0f) * M3D_PI * r * r * r;
			float i = 0.4f * mass * r * r;
			inertia = m3d_v3(i, i, i);
			break;
		}
		case M3D_SHAPE_BOX:
		{
			m3d_vec3 e = def->halfExtents;
			mass = density * 8.0f * e.x * e.y * e.z;
			inertia = m3d_v3(mass / 3.0f * (e.y * e.y + e.z * e.z), mass / 3.0f * (e.x * e.x + e.z * e.z),
							 mass / 3.0f * (e.x * e.x + e.y * e.y));
			break;
		}
		case M3D_SHAPE_CAPSULE:
		{
			float r = def->radius;
			float h = def->halfHeight;
			float mc = density * M3D_PI * r * r * (2.0f * h);
			float ms = density * (4.0f / 3.0f) * M3D_PI * r * r * r;
			mass = mc + ms;
			float iy = 0.5f * mc * r * r + 0.4f * ms * r * r;
			float d = h + 0.375f * r;
			float ixz = mc * (0.25f * r * r + h * h / 3.0f) + (83.0f / 320.0f) * ms * r * r + ms * d * d;
			inertia = m3d_v3(ixz, iy, ixz);
			break;
		}
	}

	body.mass = mass;
	if (body.type == M3D_BODY_DYNAMIC && mass > 0.0f)
	{
		body.invMass = 1.0f / mass;
		body.invInertiaLocal = m3d_v3(inertia.x > 0.0f ? 1.0f / inertia.x : 0.0f,
									  inertia.y > 0.0f ? 1.0f / inertia.y : 0.0f,
									  inertia.z > 0.0f ? 1.0f / inertia.z : 0.0f);
	}
	else
	{
		body.invMass = 0.0f;
		body.invInertiaLocal = m3d_v3_zero();
	}
}

void RemoveContactAt(m3d_world* world, int32_t index)
{
	Contact& c = world->contacts[index];
	world->contactMap.erase(c.key);
	int32_t last = (int32_t)world->contacts.size() - 1;
	if (index != last)
	{
		world->contacts[index] = world->contacts[last];
		world->contactMap[world->contacts[index].key] = index;
	}
	world->contacts.pop_back();
}

void WakeBodyInternal(m3d_world* world, int32_t index)
{
	Body& b = world->bodies[index];
	if (b.type == M3D_BODY_STATIC)
	{
		return;
	}
	if (!b.isAwake)
	{
		world->stableDirty = true;
	}
	b.isAwake = true;
	b.sleepTime = 0.0f;
}

// expanded AABB: shape AABB + margin + directional velocity extension
m3d_aabb ExpandedAABB(const Body& b, float dt)
{
	m3d_aabb aabb = ComputeShapeAABB(b.shape, BodyTransform(b));
	m3d_vec3 margin = m3d_v3(kAABBMargin, kAABBMargin, kAABBMargin);
	aabb.lowerBound = m3d_sub(aabb.lowerBound, margin);
	aabb.upperBound = m3d_add(aabb.upperBound, margin);
	m3d_vec3 d = m3d_scale(dt, b.linearVelocity);
	if (d.x < 0.0f) aabb.lowerBound.x += d.x; else aabb.upperBound.x += d.x;
	if (d.y < 0.0f) aabb.lowerBound.y += d.y; else aabb.upperBound.y += d.y;
	if (d.z < 0.0f) aabb.lowerBound.z += d.z; else aabb.upperBound.z += d.z;
	return aabb;
}

void RefreshBroadphase(m3d_world* world, float dt)
{
	const int32_t n = (int32_t)world->bodies.size();
	world->expandedAABBs.resize(n);
	world->isActiveBody.resize(n);
	world->isStableBody.resize(n);
	world->inStableTier.resize(n, 0);
	world->hasFatAABB.resize(n, 0);
	world->movedFlag.resize(n, 0);

	// stable rebuilds are coalesced: fresh sleepers keep riding the active
	// tier for a few steps so sleep/wake churn does not rebuild the stable
	// tier every step
	bool rebuildStable = world->stableForce || (world->stableDirty && world->stepsSinceStableBuild >= 4);
	world->pool->ParallelFor(n, 512, [&](int32_t i) {
		const Body& b = world->bodies[i];
		bool stable = b.inUse && !b.isAwake;
		bool active = b.inUse && (b.isAwake || !world->inStableTier[i]);
		world->isActiveBody[i] = active ? 1 : 0;
		world->isStableBody[i] = stable ? 1 : 0;
		world->movedFlag[i] = 0;
		if (active)
		{
			// fat AABB persists until the ONE-STEP-AHEAD prediction breaches it;
			// only then is the body re-paired (Box2D-style moved-proxy logic).
			// Testing the prediction, not the tight box, guarantees the pair
			// exists at least one step before a fast body can reach an obstacle.
			m3d_aabb pred = ComputeShapeAABB(b.shape, BodyTransform(b));
			m3d_vec3 step = m3d_scale(dt, b.linearVelocity);
			if (step.x < 0.0f) pred.lowerBound.x += step.x; else pred.upperBound.x += step.x;
			if (step.y < 0.0f) pred.lowerBound.y += step.y; else pred.upperBound.y += step.y;
			if (step.z < 0.0f) pred.lowerBound.z += step.z; else pred.upperBound.z += step.z;
			if (!world->hasFatAABB[i] || !m3d_aabb_contains(world->expandedAABBs[i], pred))
			{
				m3d_vec3 margin = m3d_v3(kAABBMargin, kAABBMargin, kAABBMargin);
				m3d_aabb fat = { m3d_sub(pred.lowerBound, margin), m3d_add(pred.upperBound, margin) };
				m3d_vec3 slack = m3d_scale(2.0f * dt, b.linearVelocity);
				if (slack.x < 0.0f) fat.lowerBound.x += slack.x; else fat.upperBound.x += slack.x;
				if (slack.y < 0.0f) fat.lowerBound.y += slack.y; else fat.upperBound.y += slack.y;
				if (slack.z < 0.0f) fat.lowerBound.z += slack.z; else fat.upperBound.z += slack.z;
				world->expandedAABBs[i] = fat;
				world->hasFatAABB[i] = 1;
				world->movedFlag[i] = 1;
			}
		}
		else if (stable && rebuildStable && !world->hasFatAABB[i])
		{
			world->expandedAABBs[i] = ExpandedAABB(b, dt);
			world->hasFatAABB[i] = 1;
		}
	});
	world->grid.BuildActive(world->def.gridCellSize, world->expandedAABBs, world->isActiveBody);
	if (rebuildStable)
	{
		world->grid.BuildStable(world->def.gridCellSize, world->expandedAABBs, world->isStableBody);
		for (int32_t i = 0; i < n; ++i)
		{
			world->inStableTier[i] = world->isStableBody[i] ? 1 : 0;
		}
		world->stableDirty = false;
		world->stableForce = false;
		world->stepsSinceStableBuild = 0;
	}
	else
	{
		++world->stepsSinceStableBuild;
	}
	world->gridDirty = false;
}

} // namespace

// ---------------------------------------------------------------------------
// Thread pool
// ---------------------------------------------------------------------------

namespace m3d
{

ThreadPool::ThreadPool(int32_t workerCount)
{
	int32_t extra = workerCount - 1;
	if (extra < 0)
	{
		extra = 0;
	}
	for (int32_t i = 0; i < extra; ++i)
	{
		m_threads.emplace_back(&ThreadPool::WorkerLoop, this);
	}
}

static inline void CpuRelax()
{
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
	__builtin_ia32_pause();
#else
	std::this_thread::yield();
#endif
}

ThreadPool::~ThreadPool()
{
	m_shutdown.store(true, std::memory_order_release);
	{
		std::lock_guard<std::mutex> lock(m_mutex);
	}
	m_cvWork.notify_all();
	for (std::thread& t : m_threads)
	{
		t.join();
	}
}

void ThreadPool::RunChunks()
{
	// Safe against stale wakeups: a chunk is only claimed via m_next, and a
	// finished job has m_next >= m_count, so the (possibly dangling) fn
	// pointer of an old job is never invoked.
	const std::function<void(int32_t)>* fn = m_fn;
	int32_t count = m_count;
	int32_t grain = m_grain;
	for (;;)
	{
		int32_t begin = m_next.fetch_add(grain, std::memory_order_relaxed);
		if (begin >= count)
		{
			break;
		}
		int32_t end = begin + grain < count ? begin + grain : count;
		for (int32_t i = begin; i < end; ++i)
		{
			(*fn)(i);
		}
	}
}

void ThreadPool::WorkerLoop()
{
	uint64_t lastJob = 0;
	for (;;)
	{
		if (m_shutdown.load(std::memory_order_relaxed))
		{
			return;
		}

		// job discovery is a cheap atomic spin; ENTRY goes through the mutex
		// so a straggler can never act on a half-published or retired job
		uint64_t id = m_jobId.load(std::memory_order_acquire);
		if (id != lastJob)
		{
			lastJob = id;
			bool entered = false;
			{
				std::lock_guard<std::mutex> lock(m_mutex);
				if (m_jobOpen && m_jobId.load(std::memory_order_relaxed) == id)
				{
					m_activeWorkers.fetch_add(1, std::memory_order_acq_rel);
					entered = true;
				}
			}
			if (entered)
			{
				RunChunks();
				m_activeWorkers.fetch_sub(1, std::memory_order_acq_rel);
			}
			continue;
		}

		// stay hot briefly: solver stages arrive microseconds apart
		bool sawWork = false;
		for (int32_t spin = 0; spin < 2000; ++spin)
		{
			if (m_jobId.load(std::memory_order_acquire) != lastJob || m_shutdown.load(std::memory_order_relaxed))
			{
				sawWork = true;
				break;
			}
			CpuRelax();
		}
		if (sawWork)
		{
			continue;
		}

		std::unique_lock<std::mutex> lock(m_mutex);
		m_sleepers.fetch_add(1, std::memory_order_acq_rel);
		m_cvWork.wait(lock, [&] {
			return m_shutdown.load(std::memory_order_relaxed) || m_jobId.load(std::memory_order_acquire) != lastJob;
		});
		m_sleepers.fetch_sub(1, std::memory_order_acq_rel);
	}
}

void ThreadPool::ParallelFor(int32_t count, int32_t grain, const std::function<void(int32_t)>& fn)
{
	if (count <= 0)
	{
		return;
	}
	if (m_threads.empty() || count <= grain)
	{
		for (int32_t i = 0; i < count; ++i)
		{
			fn(i);
		}
		return;
	}

	// publish under the mutex; discovery is the lock-free jobId bump
	{
		std::lock_guard<std::mutex> lock(m_mutex);
		m_fn = &fn;
		m_count = count;
		m_grain = grain > 0 ? grain : 1;
		m_next.store(0, std::memory_order_relaxed);
		m_jobOpen = true;
		m_jobId.fetch_add(1, std::memory_order_release);
	}
	if (m_sleepers.load(std::memory_order_acquire) > 0)
	{
		m_cvWork.notify_all();
	}

	RunChunks();

	// retire the job: after this no worker can enter it any more...
	{
		std::lock_guard<std::mutex> lock(m_mutex);
		m_jobOpen = false;
	}
	// ...then drain those already inside (they exit fast: m_next is spent)
	while (m_activeWorkers.load(std::memory_order_acquire) != 0)
	{
		CpuRelax();
	}
}

} // namespace m3d

// ---------------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------------

m3d_world_def m3d_world_def_default(void)
{
	m3d_world_def def;
	def.gravity = m3d_v3(0.0f, -9.81f, 0.0f);
	def.workerCount = 1;
	def.positionIterations = 2;
	def.gridCellSize = 1.5f;
	def.enableSleep = true;
	return def;
}

m3d_body_def m3d_body_def_default(void)
{
	m3d_body_def def;
	def.type = M3D_BODY_STATIC;
	def.position = m3d_v3_zero();
	def.orientation = m3d_quat_identity();
	def.linearVelocity = m3d_v3_zero();
	def.angularVelocity = m3d_v3_zero();
	def.linearDamping = 0.0f;
	def.angularDamping = 0.05f;
	def.gravityScale = 1.0f;
	def.enableSleep = true;
	def.isAwake = true;
	def.userData = 0;
	return def;
}

m3d_shape_def m3d_shape_def_default(void)
{
	m3d_shape_def def;
	def.type = M3D_SHAPE_SPHERE;
	def.radius = 0.5f;
	def.halfHeight = 0.5f;
	def.halfExtents = m3d_v3(0.5f, 0.5f, 0.5f);
	def.density = 1000.0f;
	def.friction = 0.6f;
	def.restitution = 0.0f;
	return def;
}

// ---------------------------------------------------------------------------
// World lifecycle
// ---------------------------------------------------------------------------

m3d_world* m3d_world_create(const m3d_world_def* def)
{
	m3d_world* world = new m3d_world();
	world->def = def != nullptr ? *def : m3d_world_def_default();
	if (world->def.workerCount < 1)
	{
		world->def.workerCount = 1;
	}
	if (world->def.gridCellSize <= 0.1f)
	{
		world->def.gridCellSize = 2.5f;
	}
	world->pool = new ThreadPool(world->def.workerCount);
	memset(&world->profile, 0, sizeof(world->profile));
	return world;
}

void m3d_world_destroy(m3d_world* world)
{
	delete world->pool;
	delete world;
}

// ---------------------------------------------------------------------------
// Bodies
// ---------------------------------------------------------------------------

m3d_body_id m3d_body_create(m3d_world* world, const m3d_body_def* bodyDef, const m3d_shape_def* shapeDef)
{
	int32_t index;
	if (!world->freeBodies.empty())
	{
		index = world->freeBodies.back();
		world->freeBodies.pop_back();
	}
	else
	{
		index = (int32_t)world->bodies.size();
		world->bodies.emplace_back();
		world->bodies[index].generation = 1;
	}

	Body& b = world->bodies[index];
	uint32_t gen = b.generation;
	b = Body();
	b.generation = gen;
	b.inUse = true;

	b.type = bodyDef->type;
	b.position = bodyDef->position;
	b.rotation = m3d_quat_normalize(bodyDef->orientation);
	b.linearVelocity = bodyDef->linearVelocity;
	b.angularVelocity = bodyDef->angularVelocity;
	b.linearDamping = bodyDef->linearDamping;
	b.angularDamping = bodyDef->angularDamping;
	b.gravityScale = bodyDef->gravityScale;
	b.enableSleep = bodyDef->enableSleep && world->def.enableSleep;
	b.isAwake = bodyDef->type == M3D_BODY_STATIC ? false : bodyDef->isAwake;
	b.sleepTime = 0.0f;
	b.userData = bodyDef->userData;

	b.shape.type = shapeDef->type;
	b.shape.radius = shapeDef->radius;
	b.shape.halfHeight = shapeDef->halfHeight;
	b.shape.halfExtents = shapeDef->halfExtents;
	b.friction = shapeDef->friction;
	b.restitution = shapeDef->restitution;

	ComputeMassProperties(b, shapeDef);
	UpdateInvInertiaWorld(b);

	if (index < (int32_t)world->hasFatAABB.size())
	{
		world->hasFatAABB[index] = 0; // recycled slot: stale fat AABB
	}
	world->gridDirty = true;
	world->stableDirty = true;
	return MakeBodyId(world, index);
}

void m3d_body_destroy(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	if (b == nullptr)
	{
		return;
	}
	int32_t index = id.index;

	for (int32_t i = (int32_t)world->contacts.size() - 1; i >= 0; --i)
	{
		Contact& c = world->contacts[i];
		if (c.bodyA == index || c.bodyB == index)
		{
			int32_t other = c.bodyA == index ? c.bodyB : c.bodyA;
			if (c.touching)
			{
				WakeBodyInternal(world, other);
				world->endEvents.push_back({ MakeBodyId(world, c.bodyA), MakeBodyId(world, c.bodyB) });
			}
			RemoveContactAt(world, i);
		}
	}

	for (int32_t i = 0; i < (int32_t)world->joints.size(); ++i)
	{
		Joint& j = world->joints[i];
		if (j.inUse && (j.bodyA == index || j.bodyB == index))
		{
			WakeBodyInternal(world, j.bodyA == index ? j.bodyB : j.bodyA);
			j.inUse = false;
			j.generation++;
			world->freeJoints.push_back(i);
		}
	}

	b->inUse = false;
	b->generation++;
	world->freeBodies.push_back(index);
	world->gridDirty = true;
	world->stableDirty = true;
	world->stableForce = true; // stale stable entries must not dangle
}

bool m3d_body_is_valid(m3d_world* world, m3d_body_id id)
{
	return GetBodyChecked(world, id) != nullptr;
}

m3d_vec3 m3d_body_position(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	return b != nullptr ? b->position : m3d_v3_zero();
}

m3d_quat m3d_body_rotation(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	return b != nullptr ? b->rotation : m3d_quat_identity();
}

m3d_transform m3d_body_transform(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	if (b == nullptr)
	{
		return { m3d_v3_zero(), m3d_quat_identity() };
	}
	return BodyTransform(*b);
}

void m3d_body_set_transform(m3d_world* world, m3d_body_id id, m3d_vec3 position, m3d_quat rotation)
{
	Body* b = GetBodyChecked(world, id);
	if (b == nullptr)
	{
		return;
	}
	b->position = position;
	b->rotation = m3d_quat_normalize(rotation);
	UpdateInvInertiaWorld(*b);
	if (id.index < (int32_t)world->hasFatAABB.size())
	{
		world->hasFatAABB[id.index] = 0; // teleport: fat AABB no longer valid
	}
	world->gridDirty = true;
	world->stableDirty = true;
	WakeBodyInternal(world, id.index);
}

m3d_vec3 m3d_body_linear_velocity(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	return b != nullptr ? b->linearVelocity : m3d_v3_zero();
}

m3d_vec3 m3d_body_angular_velocity(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	return b != nullptr ? b->angularVelocity : m3d_v3_zero();
}

void m3d_body_set_linear_velocity(m3d_world* world, m3d_body_id id, m3d_vec3 v)
{
	Body* b = GetBodyChecked(world, id);
	if (b == nullptr || b->type == M3D_BODY_STATIC)
	{
		return;
	}
	b->linearVelocity = v;
	WakeBodyInternal(world, id.index);
}

void m3d_body_set_angular_velocity(m3d_world* world, m3d_body_id id, m3d_vec3 w)
{
	Body* b = GetBodyChecked(world, id);
	if (b == nullptr || b->type == M3D_BODY_STATIC)
	{
		return;
	}
	b->angularVelocity = w;
	WakeBodyInternal(world, id.index);
}

void m3d_body_apply_impulse(m3d_world* world, m3d_body_id id, m3d_vec3 impulse, m3d_vec3 worldPoint)
{
	Body* b = GetBodyChecked(world, id);
	if (b == nullptr || b->type != M3D_BODY_DYNAMIC)
	{
		return;
	}
	UpdateInvInertiaWorld(*b);
	b->linearVelocity = m3d_mul_add(b->linearVelocity, b->invMass, impulse);
	m3d_vec3 r = m3d_sub(worldPoint, b->position);
	b->angularVelocity = m3d_add(b->angularVelocity, m3d_mat3_mulv(&b->invInertiaWorld, m3d_cross(r, impulse)));
	WakeBodyInternal(world, id.index);
}

float m3d_body_mass(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	return b != nullptr ? b->mass : 0.0f;
}

bool m3d_body_is_awake(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	return b != nullptr && b->isAwake;
}

void m3d_body_wake(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	if (b != nullptr)
	{
		WakeBodyInternal(world, id.index);
	}
}

uint64_t m3d_body_user_data(m3d_world* world, m3d_body_id id)
{
	Body* b = GetBodyChecked(world, id);
	return b != nullptr ? b->userData : 0;
}

// ---------------------------------------------------------------------------
// Joints
// ---------------------------------------------------------------------------

static m3d_joint_id AllocateJoint(m3d_world* world, Joint** out)
{
	int32_t index;
	if (!world->freeJoints.empty())
	{
		index = world->freeJoints.back();
		world->freeJoints.pop_back();
	}
	else
	{
		index = (int32_t)world->joints.size();
		world->joints.emplace_back();
		world->joints[index].generation = 1;
	}
	Joint& j = world->joints[index];
	uint32_t gen = j.generation;
	j = Joint();
	j.generation = gen;
	j.inUse = true;
	*out = &j;
	return { index, gen };
}

m3d_joint_id m3d_joint_create_distance(m3d_world* world, const m3d_distance_joint_def* def)
{
	Body* a = GetBodyChecked(world, def->bodyA);
	Body* b = GetBodyChecked(world, def->bodyB);
	if (a == nullptr || b == nullptr)
	{
		return { -1, 0 };
	}
	Joint* j;
	m3d_joint_id id = AllocateJoint(world, &j);
	j->type = JointType::Distance;
	j->bodyA = def->bodyA.index;
	j->bodyB = def->bodyB.index;
	j->localAnchorA = def->localAnchorA;
	j->localAnchorB = def->localAnchorB;
	j->hertz = def->hertz;
	j->dampingRatio = def->dampingRatio;
	j->lambda = 0.0f;
	if (def->length > 0.0f)
	{
		j->length = def->length;
	}
	else
	{
		m3d_vec3 pA = m3d_transform_point(BodyTransform(*a), def->localAnchorA);
		m3d_vec3 pB = m3d_transform_point(BodyTransform(*b), def->localAnchorB);
		j->length = m3d_length(m3d_sub(pB, pA));
	}
	WakeBodyInternal(world, j->bodyA);
	WakeBodyInternal(world, j->bodyB);
	return id;
}

m3d_joint_id m3d_joint_create_ball(m3d_world* world, const m3d_ball_joint_def* def)
{
	Body* a = GetBodyChecked(world, def->bodyA);
	Body* b = GetBodyChecked(world, def->bodyB);
	if (a == nullptr || b == nullptr)
	{
		return { -1, 0 };
	}
	Joint* j;
	m3d_joint_id id = AllocateJoint(world, &j);
	j->type = JointType::Ball;
	j->bodyA = def->bodyA.index;
	j->bodyB = def->bodyB.index;
	j->localAnchorA = def->localAnchorA;
	j->localAnchorB = def->localAnchorB;
	j->lambda = 0.0f;
	WakeBodyInternal(world, j->bodyA);
	WakeBodyInternal(world, j->bodyB);
	return id;
}

void m3d_joint_destroy(m3d_world* world, m3d_joint_id id)
{
	if (id.index < 0 || id.index >= (int32_t)world->joints.size())
	{
		return;
	}
	Joint& j = world->joints[id.index];
	if (!j.inUse || j.generation != id.generation)
	{
		return;
	}
	WakeBodyInternal(world, j.bodyA);
	WakeBodyInternal(world, j.bodyB);
	j.inUse = false;
	j.generation++;
	world->freeJoints.push_back(id.index);
}

// ---------------------------------------------------------------------------
// Ray casting: 3D DDA over grid cells + coarse list
// ---------------------------------------------------------------------------

m3d_ray_result m3d_world_ray_cast(m3d_world* world, m3d_vec3 origin, m3d_vec3 translation)
{
	m3d_ray_result result;
	result.body = m3d_body_id_null();
	result.point = m3d_add(origin, translation);
	result.normal = m3d_v3_zero();
	result.fraction = 2.0f;
	result.hit = false;

	if (world->gridDirty || world->stableDirty)
	{
		RefreshBroadphase(world, 0.0f);
	}

	world->rayStamps.resize(world->bodies.size(), 0);
	uint32_t stamp = ++world->rayStampCounter;

	float maxFraction = 1.0f;
	auto testBody = [&](int32_t bodyIndex) {
		if (world->rayStamps[bodyIndex] == stamp)
		{
			return;
		}
		world->rayStamps[bodyIndex] = stamp;
		Body& b = world->bodies[bodyIndex];
		if (!b.inUse)
		{
			return;
		}
		float fraction;
		m3d_vec3 normal;
		if (RayCastShape(b.shape, BodyTransform(b), origin, translation, &fraction, &normal) && fraction < maxFraction)
		{
			maxFraction = fraction;
			result.hit = true;
			result.body = MakeBodyId(world, bodyIndex);
			result.fraction = fraction;
			result.point = m3d_mul_add(origin, fraction, translation);
			result.normal = normal;
		}
	};

	// coarse lists first (large statics such as the ground)
	for (int32_t bodyIndex : world->grid.ActiveOversized())
	{
		testBody(bodyIndex);
	}
	for (int32_t bodyIndex : world->grid.StableOversized())
	{
		testBody(bodyIndex);
	}

	// walk grid cells along the segment (Amanatides & Woo)
	const auto& activeEntries = world->grid.ActiveEntries();
	const auto& stableEntries = world->grid.StableEntries();
	if (!activeEntries.empty() || !stableEntries.empty())
	{
		float cs = world->grid.CellSize();
		float inv = 1.0f / cs;
		int64_t ix = (int64_t)floorf(origin.x * inv);
		int64_t iy = (int64_t)floorf(origin.y * inv);
		int64_t iz = (int64_t)floorf(origin.z * inv);

		const float* d = (const float*)&translation;
		const float* o = (const float*)&origin;
		int64_t cell[3] = { ix, iy, iz };
		int64_t stepDir[3];
		float tMax[3], tDelta[3];
		for (int k = 0; k < 3; ++k)
		{
			if (fabsf(d[k]) < 1e-12f)
			{
				stepDir[k] = 0;
				tMax[k] = FLT_MAX;
				tDelta[k] = FLT_MAX;
			}
			else
			{
				stepDir[k] = d[k] > 0.0f ? 1 : -1;
				float boundary = (float)(cell[k] + (stepDir[k] > 0 ? 1 : 0)) * cs;
				tMax[k] = (boundary - o[k]) / d[k];
				tDelta[k] = cs / fabsf(d[k]);
			}
		}

		auto visitCell = [&](int64_t x, int64_t y, int64_t z) {
			// binary search both sorted tiers for this cell's key
			m3d_vec3 probe = m3d_v3(((float)x + 0.5f) * cs, ((float)y + 0.5f) * cs, ((float)z + 0.5f) * cs);
			uint64_t key = world->grid.CellKeyOf(probe);
			for (const auto* entries : { &activeEntries, &stableEntries })
			{
				auto it = std::lower_bound(entries->begin(), entries->end(), key,
										   [](const SpatialGrid::Entry& e, uint64_t k) { return e.key < k; });
				while (it != entries->end() && it->key == key)
				{
					testBody(it->body);
					++it;
				}
			}
		};

		float t = 0.0f;
		for (int guard = 0; guard < 4096 && t <= maxFraction; ++guard)
		{
			visitCell(cell[0], cell[1], cell[2]);
			int axis = 0;
			if (tMax[1] < tMax[axis]) axis = 1;
			if (tMax[2] < tMax[axis]) axis = 2;
			t = tMax[axis];
			if (t > 1.0f || stepDir[axis] == 0)
			{
				break;
			}
			cell[axis] += stepDir[axis];
			tMax[axis] += tDelta[axis];
		}
	}

	return result;
}

m3d_contact_events m3d_world_contact_events(m3d_world* world)
{
	m3d_contact_events events;
	events.beginEvents = world->beginEvents.data();
	events.beginCount = (int32_t)world->beginEvents.size();
	events.endEvents = world->endEvents.data();
	events.endCount = (int32_t)world->endEvents.size();
	return events;
}

m3d_profile m3d_world_profile(m3d_world* world)
{
	return world->profile;
}

// ---------------------------------------------------------------------------
// Step
// ---------------------------------------------------------------------------

void m3d_world_step(m3d_world* world, float dt, int substepCount)
{
	if (dt <= 0.0f)
	{
		return;
	}
	if (substepCount < 1)
	{
		substepCount = 1;
	}

	auto tStart = Clock::now();
	world->beginEvents.clear();
	world->endEvents.clear();

	const int32_t bodyCapacity = (int32_t)world->bodies.size();

	// --- broadphase: rebuild the hash grid over expanded AABBs ---------------
	RefreshBroadphase(world, dt);

	// pairs connected by a joint do not collide
	std::unordered_set<uint64_t> jointedPairs;
	for (const Joint& j : world->joints)
	{
		if (j.inUse)
		{
			int32_t lo = j.bodyA < j.bodyB ? j.bodyA : j.bodyB;
			int32_t hi = j.bodyA < j.bodyB ? j.bodyB : j.bodyA;
			jointedPairs.insert(ContactKey(lo, hi));
		}
	}
	const bool haveJoints = !jointedPairs.empty();

	// --- candidate pairs from the grid ---------------------------------------
	world->grid.QueryPairs(world->expandedAABBs, [&](int32_t a, int32_t b) {
		if (a == b)
		{
			return; // body present in both tiers between stable rebuilds
		}
		if (!world->movedFlag[a] && !world->movedFlag[b])
		{
			return; // both inside their fat AABBs: the pair set cannot have changed
		}
		int32_t lo = a < b ? a : b;
		int32_t hi = a < b ? b : a;
		Body& lb = world->bodies[lo];
		Body& hb = world->bodies[hi];
		if (!lb.inUse || !hb.inUse)
		{
			return;
		}
		if (lb.type != M3D_BODY_DYNAMIC && hb.type != M3D_BODY_DYNAMIC)
		{
			return;
		}
		if (!lb.isAwake && !hb.isAwake)
		{
			return;
		}
		if (!m3d_aabb_overlap(world->expandedAABBs[lo], world->expandedAABBs[hi]))
		{
			return;
		}
		uint64_t key = ContactKey(lo, hi);
		if (world->contactMap.count(key) != 0 || (haveJoints && jointedPairs.count(key) != 0))
		{
			return;
		}
		Contact c;
		c.bodyA = lo;
		c.bodyB = hi;
		c.key = key;
		c.manifold.pointCount = 0;
		c.manifold.normal = m3d_v3(0.0f, 1.0f, 0.0f);
		c.friction = sqrtf(lb.friction * hb.friction);
		c.restitution = lb.restitution > hb.restitution ? lb.restitution : hb.restitution;
		c.touching = false;
		c.hasCache = false;
		world->contactMap[key] = (int32_t)world->contacts.size();
		world->contacts.push_back(c);
	});

	// --- destroy contacts whose expanded AABBs separated ----------------------
	// parallel read-only scan, then serial removal in descending index order
	{
		const int32_t scanCount = (int32_t)world->contacts.size();
		std::vector<uint8_t> removeFlag(scanCount, 0);
		world->pool->ParallelFor(scanCount, 512, [&](int32_t i) {
			Contact& c = world->contacts[i];
			Body& a = world->bodies[c.bodyA];
			Body& b = world->bodies[c.bodyB];
			if (!a.isAwake && !b.isAwake)
			{
				return;
			}
			if (!m3d_aabb_overlap(world->expandedAABBs[c.bodyA], world->expandedAABBs[c.bodyB]))
			{
				removeFlag[i] = 1;
			}
		});
		for (int32_t i = scanCount - 1; i >= 0; --i)
		{
			if (!removeFlag[i])
			{
				continue;
			}
			Contact& c = world->contacts[i];
			if (c.touching)
			{
				world->endEvents.push_back({ MakeBodyId(world, c.bodyA), MakeBodyId(world, c.bodyB) });
			}
			RemoveContactAt(world, i);
		}
	}
	auto tBroad = Clock::now();

	// --- narrowphase (parallel only when there is real collision work) --------
	// cached-manifold steps cost ~30ns per contact: waking workers for that is
	// slower than doing it inline. Last step's regeneration count is the gate.
	const int32_t contactCount = (int32_t)world->contacts.size();
	std::vector<int8_t> transition(contactCount, 0);
	std::atomic<int32_t> regenCount{ 0 };

	bool parallelNp = world->pool->WorkerCount() > 1 && world->lastNpRegenCount > 256;
	int32_t npGrain = contactCount / (4 * world->pool->WorkerCount()) + 1;
	npGrain = npGrain < 64 ? 64 : npGrain;
	if (!parallelNp)
	{
		npGrain = contactCount + 1; // forces the inline path
	}
	world->pool->ParallelFor(contactCount, npGrain, [&](int32_t i) {
		Contact& c = world->contacts[i];
		Body& a = world->bodies[c.bodyA];
		Body& b = world->bodies[c.bodyB];
		if (!a.isAwake && !b.isAwake)
		{
			return;
		}

		// reuse the cached manifold while neither body has moved
		if (!c.hasCache || PoseChanged(a, c.cachePosA, c.cacheRotA) || PoseChanged(b, c.cachePosB, c.cacheRotB))
		{
			regenCount.fetch_add(1, std::memory_order_relaxed);
			Manifold oldManifold = c.manifold;
			// speculative reach scales with the relative approach this step
			float relSpeed = m3d_length(m3d_sub(a.linearVelocity, b.linearVelocity));
			float margin = kSpeculativeBase + dt * relSpeed;
			CollideShapes(a.shape, BodyTransform(a), b.shape, BodyTransform(b), margin, &c.manifold);
			InheritFrictionAnchors(oldManifold, c.manifold, a, b);
			c.hasCache = true;
			c.cachePosA = a.position;
			c.cacheRotA = a.rotation;
			c.cachePosB = b.position;
			c.cacheRotB = b.rotation;
			for (int32_t p = 0; p < c.manifold.pointCount; ++p)
			{
				ManifoldPoint& mp = c.manifold.points[p];
				m3d_vec3 rA = m3d_rotate(a.rotation, mp.localAnchorA);
				m3d_vec3 rB = m3d_rotate(b.rotation, mp.localAnchorB);
				float ref = m3d_dot(m3d_sub(m3d_add(b.position, rB), m3d_add(a.position, rA)), c.manifold.normal);
				mp.separationOffset = mp.separation - ref;
			}
		}

		bool wasTouching = c.touching;
		c.touching = c.manifold.pointCount > 0;
		if (c.touching && !wasTouching)
		{
			transition[i] = 1;
		}
		else if (!c.touching && wasTouching)
		{
			transition[i] = -1;
		}
	});

	world->lastNpRegenCount = regenCount.load(std::memory_order_relaxed);

	for (int32_t i = 0; i < contactCount; ++i)
	{
		if (transition[i] == 0)
		{
			continue;
		}
		Contact& c = world->contacts[i];
		if (transition[i] > 0)
		{
			world->beginEvents.push_back({ MakeBodyId(world, c.bodyA), MakeBodyId(world, c.bodyB) });
			Body& a = world->bodies[c.bodyA];
			Body& b = world->bodies[c.bodyB];
			if (a.isAwake && !b.isAwake)
			{
				WakeBodyInternal(world, c.bodyB);
			}
			else if (b.isAwake && !a.isAwake)
			{
				WakeBodyInternal(world, c.bodyA);
			}
		}
		else
		{
			world->endEvents.push_back({ MakeBodyId(world, c.bodyA), MakeBodyId(world, c.bodyB) });
		}
	}
	auto tNarrow = Clock::now();

	// --- islands via union-find over dynamic bodies ----------------------------
	std::vector<int32_t> parent(bodyCapacity);
	for (int32_t i = 0; i < bodyCapacity; ++i)
	{
		parent[i] = i;
	}
	auto find = [&](int32_t x) {
		while (parent[x] != x)
		{
			parent[x] = parent[parent[x]];
			x = parent[x];
		}
		return x;
	};
	auto unite = [&](int32_t x, int32_t y) {
		int32_t rx = find(x);
		int32_t ry = find(y);
		if (rx != ry)
		{
			parent[rx] = ry;
		}
	};

	for (const Contact& c : world->contacts)
	{
		if (!c.touching)
		{
			continue;
		}
		if (world->bodies[c.bodyA].type == M3D_BODY_DYNAMIC && world->bodies[c.bodyB].type == M3D_BODY_DYNAMIC)
		{
			unite(c.bodyA, c.bodyB);
		}
	}
	for (const Joint& j : world->joints)
	{
		if (j.inUse && world->bodies[j.bodyA].type == M3D_BODY_DYNAMIC && world->bodies[j.bodyB].type == M3D_BODY_DYNAMIC)
		{
			unite(j.bodyA, j.bodyB);
		}
	}

	std::vector<int32_t> islandIndex(bodyCapacity, -1);
	world->islandBodies.clear();
	world->islandContacts.clear();
	world->islandJoints.clear();

	std::vector<bool> islandAwake;
	for (int32_t i = 0; i < bodyCapacity; ++i)
	{
		Body& b = world->bodies[i];
		if (!b.inUse || b.type != M3D_BODY_DYNAMIC)
		{
			continue;
		}
		int32_t root = find(i);
		if (islandIndex[root] < 0)
		{
			islandIndex[root] = (int32_t)world->islandBodies.size();
			world->islandBodies.emplace_back();
			world->islandContacts.emplace_back();
			world->islandJoints.emplace_back();
			islandAwake.push_back(false);
		}
		int32_t idx = islandIndex[root];
		world->islandBodies[idx].push_back(i);
		if (b.isAwake)
		{
			islandAwake[idx] = true;
		}
	}

	const int32_t islandCount = (int32_t)world->islandBodies.size();
	for (int32_t idx = 0; idx < islandCount; ++idx)
	{
		if (!islandAwake[idx])
		{
			continue;
		}
		for (int32_t bi : world->islandBodies[idx])
		{
			if (!world->bodies[bi].isAwake)
			{
				WakeBodyInternal(world, bi);
			}
		}
	}

	int32_t touchingCount = 0;
	for (int32_t i = 0; i < (int32_t)world->contacts.size(); ++i)
	{
		Contact& c = world->contacts[i];
		if (!c.touching)
		{
			continue;
		}
		++touchingCount;
		int32_t dynBody = world->bodies[c.bodyA].type == M3D_BODY_DYNAMIC ? c.bodyA : c.bodyB;
		int32_t idx = islandIndex[find(dynBody)];
		if (idx >= 0 && islandAwake[idx])
		{
			world->islandContacts[idx].push_back(i);
		}
	}
	for (int32_t i = 0; i < (int32_t)world->joints.size(); ++i)
	{
		Joint& j = world->joints[i];
		if (!j.inUse)
		{
			continue;
		}
		int32_t dynBody = world->bodies[j.bodyA].type == M3D_BODY_DYNAMIC ? j.bodyA : j.bodyB;
		if (world->bodies[dynBody].type != M3D_BODY_DYNAMIC)
		{
			continue;
		}
		int32_t idx = islandIndex[find(dynBody)];
		if (idx >= 0 && islandAwake[idx])
		{
			world->islandJoints[idx].push_back(i);
		}
	}

	// --- solve islands in parallel (islands touch disjoint dynamic bodies) ----
	std::vector<int32_t> awakeIslands;
	for (int32_t idx = 0; idx < islandCount; ++idx)
	{
		if (islandAwake[idx])
		{
			awakeIslands.push_back(idx);
		}
	}
	std::sort(awakeIslands.begin(), awakeIslands.end(), [&](int32_t l, int32_t r) {
		return world->islandBodies[l].size() > world->islandBodies[r].size();
	});

	// big islands first: solved one at a time with graph coloring, which
	// parallelizes INSIDE the island. The colored path is chosen by island
	// size only (never by thread count), so results stay bit-identical for
	// any workerCount: within a color bodies are disjoint.
	std::vector<int32_t> smallIslands;
	smallIslands.reserve(awakeIslands.size());
	for (int32_t idx : awakeIslands)
	{
		if ((int32_t)world->islandContacts[idx].size() >= kBigIslandContacts)
		{
			SolveIslandColoredXPBD(world, world->islandBodies[idx], world->islandContacts[idx],
								   world->islandJoints[idx], dt, substepCount);
		}
		else
		{
			smallIslands.push_back(idx);
		}
	}

	// grain 1: islands are sorted biggest-first, so batching would hand all
	// the big islands to a single worker; per-island dispatch is one atomic.
	// Near-idle worlds (a few tiny islands) run inline instead.
	int32_t awakeContactTotal = 0;
	for (int32_t idx : smallIslands)
	{
		awakeContactTotal += (int32_t)world->islandContacts[idx].size();
	}
	int32_t islandGrain = awakeContactTotal < 64 ? (int32_t)smallIslands.size() + 1 : 1;
	world->pool->ParallelFor((int32_t)smallIslands.size(), islandGrain, [&](int32_t k) {
		int32_t idx = smallIslands[k];
		SolveIslandXPBD(world, world->islandBodies[idx], world->islandContacts[idx], world->islandJoints[idx], dt,
						substepCount);
	});

	// integrate kinematic bodies
	for (int32_t i = 0; i < bodyCapacity; ++i)
	{
		Body& b = world->bodies[i];
		if (!b.inUse || b.type != M3D_BODY_KINEMATIC || !b.isAwake)
		{
			continue;
		}
		b.position = m3d_mul_add(b.position, dt, b.linearVelocity);
		b.rotation = m3d_quat_integrate(b.rotation, b.angularVelocity, dt);
	}

	// --- sleeping ---------------------------------------------------------------
	int32_t awakeBodyCount = 0;
	for (int32_t idx : awakeIslands)
	{
		float minSleepTime = FLT_MAX;
		bool canSleep = true;
		for (int32_t bi : world->islandBodies[idx])
		{
			Body& b = world->bodies[bi];
			if (!b.enableSleep)
			{
				canSleep = false;
				break;
			}
			if (m3d_length_sq(b.linearVelocity) > kSleepLinearVelSq ||
				m3d_length_sq(b.angularVelocity) > kSleepAngularVelSq)
			{
				b.sleepTime = 0.0f;
			}
			else
			{
				b.sleepTime += dt;
			}
			minSleepTime = fminf(minSleepTime, b.sleepTime);
		}
		if (canSleep && minSleepTime >= kTimeToSleep)
		{
			for (int32_t bi : world->islandBodies[idx])
			{
				Body& b = world->bodies[bi];
				b.isAwake = false;
				b.linearVelocity = m3d_v3_zero();
				b.angularVelocity = m3d_v3_zero();
			}
			world->stableDirty = true;
		}
	}

	int32_t bodyCount = 0;
	for (int32_t i = 0; i < bodyCapacity; ++i)
	{
		Body& b = world->bodies[i];
		if (!b.inUse)
		{
			continue;
		}
		++bodyCount;
		if (b.isAwake)
		{
			++awakeBodyCount;
		}
	}

	auto tEnd = Clock::now();
	world->profile.broadphaseMs = ToMs(tStart, tBroad);
	world->profile.narrowphaseMs = ToMs(tBroad, tNarrow);
	world->profile.solverMs = ToMs(tNarrow, tEnd);
	world->profile.stepMs = ToMs(tStart, tEnd);
	world->profile.bodyCount = bodyCount;
	world->profile.awakeBodyCount = awakeBodyCount;
	world->profile.contactCount = (int32_t)world->contacts.size();
	world->profile.touchingContactCount = touchingCount;
	world->profile.islandCount = islandCount;
}
