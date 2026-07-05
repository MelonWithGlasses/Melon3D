// Melon3D 1.0 - a compact 3D rigid body physics engine for games, built on
// position based dynamics (XPBD).
//
//  - Shapes: sphere, capsule, box (OBB)
//  - Broadphase: two-tier sorted spatial hash grid (sleeping/static bodies
//    live in a stable tier), moved-only pair updates via persistent fat
//    AABBs with one-step velocity prediction (anti-tunneling)
//  - Narrowphase: pose-delta manifold cache; speculative margins
//  - Solver: XPBD substepping (Mueller et al. 2020) with persistent
//    stiction friction anchors, per-substep frozen masses and angular
//    Jacobians. Zero rest penetration by construction
//  - Islands with whole-island sleeping; graph-colored solve with 8-wide
//    SIMD contact packets for big islands
//  - Joints: distance (rigid or spring via compliance), ball (spherical)
//  - Built-in spin fork-join thread pool; results are bit-identical for
//    any thread count
//  - Contact begin/end events, DDA grid ray casting, per-stage profiling
//
// SPDX-License-Identifier: MIT

#pragma once

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define M3D_VERSION_MAJOR 1
#define M3D_VERSION_MINOR 0
#define M3D_VERSION_PATCH 0

#define M3D_PI 3.14159265358979323846f

#ifdef __cplusplus
#define M3D_API extern
#else
#define M3D_API
#endif

// ---------------------------------------------------------------------------
// Math types
// ---------------------------------------------------------------------------

typedef struct m3d_vec3
{
	float x, y, z;
} m3d_vec3;

typedef struct m3d_quat
{
	float x, y, z, w;
} m3d_quat;

// Column-major 3x3 matrix
typedef struct m3d_mat3
{
	m3d_vec3 cx, cy, cz;
} m3d_mat3;

typedef struct m3d_transform
{
	m3d_vec3 p;
	m3d_quat q;
} m3d_transform;

typedef struct m3d_aabb
{
	m3d_vec3 lowerBound;
	m3d_vec3 upperBound;
} m3d_aabb;

// ---------------------------------------------------------------------------
// Inline math
// ---------------------------------------------------------------------------

static inline m3d_vec3 m3d_v3(float x, float y, float z)
{
	m3d_vec3 v = { x, y, z };
	return v;
}

static inline m3d_vec3 m3d_v3_zero(void) { return m3d_v3(0.0f, 0.0f, 0.0f); }

static inline m3d_vec3 m3d_add(m3d_vec3 a, m3d_vec3 b) { return m3d_v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline m3d_vec3 m3d_sub(m3d_vec3 a, m3d_vec3 b) { return m3d_v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline m3d_vec3 m3d_neg(m3d_vec3 a) { return m3d_v3(-a.x, -a.y, -a.z); }
static inline m3d_vec3 m3d_scale(float s, m3d_vec3 v) { return m3d_v3(s * v.x, s * v.y, s * v.z); }
static inline m3d_vec3 m3d_mul_add(m3d_vec3 a, float s, m3d_vec3 b) { return m3d_v3(a.x + s * b.x, a.y + s * b.y, a.z + s * b.z); }
static inline float m3d_dot(m3d_vec3 a, m3d_vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline m3d_vec3 m3d_cross(m3d_vec3 a, m3d_vec3 b)
{
	return m3d_v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
static inline float m3d_length_sq(m3d_vec3 v) { return m3d_dot(v, v); }
static inline float m3d_length(m3d_vec3 v) { return sqrtf(m3d_dot(v, v)); }
static inline m3d_vec3 m3d_normalize(m3d_vec3 v)
{
	float len = m3d_length(v);
	if (len < 1e-12f)
	{
		return m3d_v3_zero();
	}
	float inv = 1.0f / len;
	return m3d_scale(inv, v);
}
static inline m3d_vec3 m3d_min(m3d_vec3 a, m3d_vec3 b)
{
	return m3d_v3(a.x < b.x ? a.x : b.x, a.y < b.y ? a.y : b.y, a.z < b.z ? a.z : b.z);
}
static inline m3d_vec3 m3d_max(m3d_vec3 a, m3d_vec3 b)
{
	return m3d_v3(a.x > b.x ? a.x : b.x, a.y > b.y ? a.y : b.y, a.z > b.z ? a.z : b.z);
}
static inline m3d_vec3 m3d_abs(m3d_vec3 a) { return m3d_v3(fabsf(a.x), fabsf(a.y), fabsf(a.z)); }
static inline m3d_vec3 m3d_lerp(m3d_vec3 a, m3d_vec3 b, float t) { return m3d_add(a, m3d_scale(t, m3d_sub(b, a))); }

static inline m3d_quat m3d_quat_identity(void)
{
	m3d_quat q = { 0.0f, 0.0f, 0.0f, 1.0f };
	return q;
}

static inline m3d_quat m3d_quat_axis_angle(m3d_vec3 axis, float angle)
{
	m3d_vec3 n = m3d_normalize(axis);
	float half = 0.5f * angle;
	float s = sinf(half);
	m3d_quat q = { s * n.x, s * n.y, s * n.z, cosf(half) };
	return q;
}

static inline m3d_quat m3d_quat_mul(m3d_quat a, m3d_quat b)
{
	m3d_quat q;
	q.x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
	q.y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
	q.z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;
	q.w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
	return q;
}

static inline m3d_quat m3d_quat_normalize(m3d_quat q)
{
	float len = sqrtf(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
	if (len < 1e-12f)
	{
		return m3d_quat_identity();
	}
	float inv = 1.0f / len;
	m3d_quat r = { q.x * inv, q.y * inv, q.z * inv, q.w * inv };
	return r;
}

static inline m3d_quat m3d_quat_conj(m3d_quat q)
{
	m3d_quat r = { -q.x, -q.y, -q.z, q.w };
	return r;
}

static inline m3d_vec3 m3d_rotate(m3d_quat q, m3d_vec3 v)
{
	// v' = v + 2 * cross(q.xyz, cross(q.xyz, v) + q.w * v)
	m3d_vec3 u = m3d_v3(q.x, q.y, q.z);
	m3d_vec3 t = m3d_scale(2.0f, m3d_cross(u, v));
	return m3d_add(m3d_add(v, m3d_scale(q.w, t)), m3d_cross(u, t));
}

static inline m3d_vec3 m3d_inv_rotate(m3d_quat q, m3d_vec3 v)
{
	return m3d_rotate(m3d_quat_conj(q), v);
}

// Integrate rotation: q' = normalize(q + 0.5 * h * omega * q)
static inline m3d_quat m3d_quat_integrate(m3d_quat q, m3d_vec3 omega, float h)
{
	m3d_quat wq = { omega.x, omega.y, omega.z, 0.0f };
	m3d_quat dq = m3d_quat_mul(wq, q);
	m3d_quat r = { q.x + 0.5f * h * dq.x, q.y + 0.5f * h * dq.y, q.z + 0.5f * h * dq.z, q.w + 0.5f * h * dq.w };
	return m3d_quat_normalize(r);
}

static inline m3d_mat3 m3d_quat_to_mat3(m3d_quat q)
{
	float x = q.x, y = q.y, z = q.z, w = q.w;
	float x2 = x + x, y2 = y + y, z2 = z + z;
	float xx = x * x2, xy = x * y2, xz = x * z2;
	float yy = y * y2, yz = y * z2, zz = z * z2;
	float wx = w * x2, wy = w * y2, wz = w * z2;
	m3d_mat3 m;
	m.cx = m3d_v3(1.0f - (yy + zz), xy + wz, xz - wy);
	m.cy = m3d_v3(xy - wz, 1.0f - (xx + zz), yz + wx);
	m.cz = m3d_v3(xz + wy, yz - wx, 1.0f - (xx + yy));
	return m;
}

static inline m3d_vec3 m3d_mat3_mulv(const m3d_mat3* m, m3d_vec3 v)
{
	return m3d_add(m3d_add(m3d_scale(v.x, m->cx), m3d_scale(v.y, m->cy)), m3d_scale(v.z, m->cz));
}

// transpose(m) * v
static inline m3d_vec3 m3d_mat3_mulv_t(const m3d_mat3* m, m3d_vec3 v)
{
	return m3d_v3(m3d_dot(m->cx, v), m3d_dot(m->cy, v), m3d_dot(m->cz, v));
}

static inline m3d_vec3 m3d_transform_point(m3d_transform xf, m3d_vec3 p)
{
	return m3d_add(m3d_rotate(xf.q, p), xf.p);
}

static inline m3d_vec3 m3d_inv_transform_point(m3d_transform xf, m3d_vec3 p)
{
	return m3d_inv_rotate(xf.q, m3d_sub(p, xf.p));
}

static inline bool m3d_aabb_overlap(m3d_aabb a, m3d_aabb b)
{
	return a.lowerBound.x <= b.upperBound.x && a.lowerBound.y <= b.upperBound.y && a.lowerBound.z <= b.upperBound.z &&
		   b.lowerBound.x <= a.upperBound.x && b.lowerBound.y <= a.upperBound.y && b.lowerBound.z <= a.upperBound.z;
}

static inline bool m3d_aabb_contains(m3d_aabb a, m3d_aabb b)
{
	return a.lowerBound.x <= b.lowerBound.x && a.lowerBound.y <= b.lowerBound.y && a.lowerBound.z <= b.lowerBound.z &&
		   b.upperBound.x <= a.upperBound.x && b.upperBound.y <= a.upperBound.y && b.upperBound.z <= a.upperBound.z;
}

// ---------------------------------------------------------------------------
// Ids and definitions
// ---------------------------------------------------------------------------

typedef struct m3d_world m3d_world;

typedef struct m3d_body_id
{
	int32_t index;
	uint32_t generation;
} m3d_body_id;

typedef struct m3d_joint_id
{
	int32_t index;
	uint32_t generation;
} m3d_joint_id;

static inline m3d_body_id m3d_body_id_null(void)
{
	m3d_body_id id = { -1, 0 };
	return id;
}

typedef enum m3d_body_type
{
	M3D_BODY_STATIC = 0,
	M3D_BODY_KINEMATIC = 1,
	M3D_BODY_DYNAMIC = 2,
} m3d_body_type;

typedef enum m3d_shape_type
{
	M3D_SHAPE_SPHERE = 0,
	M3D_SHAPE_CAPSULE = 1,
	M3D_SHAPE_BOX = 2,
} m3d_shape_type;

typedef struct m3d_shape_def
{
	m3d_shape_type type;

	// sphere / capsule
	float radius;
	// capsule: half length of the inner segment (local Y axis)
	float halfHeight;
	// box
	m3d_vec3 halfExtents;

	float density;	   // kg/m^3
	float friction;	   // Coulomb friction coefficient
	float restitution; // bounciness [0,1]
} m3d_shape_def;

typedef struct m3d_body_def
{
	m3d_body_type type;
	m3d_vec3 position;
	m3d_quat orientation;
	m3d_vec3 linearVelocity;
	m3d_vec3 angularVelocity;
	float linearDamping;
	float angularDamping;
	float gravityScale;
	bool enableSleep;
	bool isAwake;
	uint64_t userData;
} m3d_body_def;

typedef struct m3d_world_def
{
	m3d_vec3 gravity;
	// Number of worker threads (including the calling thread). 1 = single threaded.
	int32_t workerCount;
	// XPBD position iterations per substep (default 2)
	int32_t positionIterations;
	// Broadphase hash grid cell size in meters (default 2.5)
	float gridCellSize;
	bool enableSleep;
} m3d_world_def;

typedef struct m3d_distance_joint_def
{
	m3d_body_id bodyA;
	m3d_body_id bodyB;
	m3d_vec3 localAnchorA;
	m3d_vec3 localAnchorB;
	float length; // <= 0 means: compute from initial pose
	float hertz;  // 0 = rigid
	float dampingRatio;
} m3d_distance_joint_def;

typedef struct m3d_ball_joint_def
{
	m3d_body_id bodyA;
	m3d_body_id bodyB;
	m3d_vec3 localAnchorA;
	m3d_vec3 localAnchorB;
} m3d_ball_joint_def;

// ---------------------------------------------------------------------------
// Events, queries, profiling
// ---------------------------------------------------------------------------

typedef struct m3d_contact_event
{
	m3d_body_id bodyA;
	m3d_body_id bodyB;
} m3d_contact_event;

typedef struct m3d_contact_events
{
	const m3d_contact_event* beginEvents;
	int32_t beginCount;
	const m3d_contact_event* endEvents;
	int32_t endCount;
} m3d_contact_events;

typedef struct m3d_ray_result
{
	m3d_body_id body;
	m3d_vec3 point;
	m3d_vec3 normal;
	float fraction; // [0,1] along the translation, or >1 if no hit
	bool hit;
} m3d_ray_result;

typedef struct m3d_profile
{
	float stepMs;
	float broadphaseMs;
	float narrowphaseMs;
	float solverMs;
	int32_t bodyCount;
	int32_t awakeBodyCount;
	int32_t contactCount;
	int32_t touchingContactCount;
	int32_t islandCount;
} m3d_profile;

// ---------------------------------------------------------------------------
// API
// ---------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#endif

M3D_API m3d_world_def m3d_world_def_default(void);
M3D_API m3d_body_def m3d_body_def_default(void);
M3D_API m3d_shape_def m3d_shape_def_default(void);

M3D_API m3d_world* m3d_world_create(const m3d_world_def* def);
M3D_API void m3d_world_destroy(m3d_world* world);

// Advance the simulation. substepCount of 4 is a good default.
M3D_API void m3d_world_step(m3d_world* world, float dt, int substepCount);

M3D_API m3d_body_id m3d_body_create(m3d_world* world, const m3d_body_def* bodyDef, const m3d_shape_def* shapeDef);
M3D_API void m3d_body_destroy(m3d_world* world, m3d_body_id id);
M3D_API bool m3d_body_is_valid(m3d_world* world, m3d_body_id id);

M3D_API m3d_vec3 m3d_body_position(m3d_world* world, m3d_body_id id);
M3D_API m3d_quat m3d_body_rotation(m3d_world* world, m3d_body_id id);
M3D_API m3d_transform m3d_body_transform(m3d_world* world, m3d_body_id id);
M3D_API void m3d_body_set_transform(m3d_world* world, m3d_body_id id, m3d_vec3 position, m3d_quat rotation);
M3D_API m3d_vec3 m3d_body_linear_velocity(m3d_world* world, m3d_body_id id);
M3D_API m3d_vec3 m3d_body_angular_velocity(m3d_world* world, m3d_body_id id);
M3D_API void m3d_body_set_linear_velocity(m3d_world* world, m3d_body_id id, m3d_vec3 v);
M3D_API void m3d_body_set_angular_velocity(m3d_world* world, m3d_body_id id, m3d_vec3 w);
M3D_API void m3d_body_apply_impulse(m3d_world* world, m3d_body_id id, m3d_vec3 impulse, m3d_vec3 worldPoint);
M3D_API float m3d_body_mass(m3d_world* world, m3d_body_id id);
M3D_API bool m3d_body_is_awake(m3d_world* world, m3d_body_id id);
M3D_API void m3d_body_wake(m3d_world* world, m3d_body_id id);
M3D_API uint64_t m3d_body_user_data(m3d_world* world, m3d_body_id id);

M3D_API m3d_joint_id m3d_joint_create_distance(m3d_world* world, const m3d_distance_joint_def* def);
M3D_API m3d_joint_id m3d_joint_create_ball(m3d_world* world, const m3d_ball_joint_def* def);
M3D_API void m3d_joint_destroy(m3d_world* world, m3d_joint_id id);

// Cast a ray from origin along translation. Returns the closest hit.
M3D_API m3d_ray_result m3d_world_ray_cast(m3d_world* world, m3d_vec3 origin, m3d_vec3 translation);

// Contact begin/end events for the last step. Valid until the next step.
M3D_API m3d_contact_events m3d_world_contact_events(m3d_world* world);

M3D_API m3d_profile m3d_world_profile(m3d_world* world);

#ifdef __cplusplus
}
#endif
