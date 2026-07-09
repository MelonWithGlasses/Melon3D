// Melon3D integration tests: falling/resting bodies, stacking stability,
// restitution, sleeping, ray casts, joints, contact events.
#include "melon3d.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int g_failCount = 0;
static int g_testCount = 0;

#define CHECK(cond, msg)                                                                                               \
	do                                                                                                                 \
	{                                                                                                                  \
		++g_testCount;                                                                                                 \
		if (!(cond))                                                                                                   \
		{                                                                                                              \
			++g_failCount;                                                                                             \
			printf("FAIL: %s (%s:%d)\n", msg, __FILE__, __LINE__);                                                     \
		}                                                                                                              \
	} while (0)

static m3d_body_id AddGround(m3d_world* world)
{
	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_STATIC;
	bd.position = m3d_v3(0.0f, -1.0f, 0.0f);
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_BOX;
	sd.halfExtents = m3d_v3(50.0f, 1.0f, 50.0f);
	return m3d_body_create(world, &bd, &sd);
}

static m3d_body_id AddBox(m3d_world* world, m3d_vec3 pos, m3d_vec3 halfExtents, float friction)
{
	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = pos;
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_BOX;
	sd.halfExtents = halfExtents;
	sd.friction = friction;
	return m3d_body_create(world, &bd, &sd);
}

// A sphere dropped on the ground must come to rest at y = radius.
static void TestSphereRest(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = m3d_v3(0.0f, 3.0f, 0.0f);
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_SPHERE;
	sd.radius = 0.5f;
	m3d_body_id ball = m3d_body_create(world, &bd, &sd);

	for (int i = 0; i < 240; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	m3d_vec3 p = m3d_body_position(world, ball);
	CHECK(fabsf(p.y - 0.5f) < 0.02f, "sphere rests at y=radius");
	CHECK(fabsf(p.x) < 0.01f && fabsf(p.z) < 0.01f, "sphere does not drift");
	m3d_vec3 v = m3d_body_linear_velocity(world, ball);
	CHECK(m3d_length(v) < 0.05f, "sphere at rest has ~zero velocity");
	m3d_world_destroy(world);
}

// A capsule lying on its side must come to rest at y = radius.
static void TestCapsuleRest(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = m3d_v3(0.0f, 2.0f, 0.0f);
	bd.orientation = m3d_quat_axis_angle(m3d_v3(0.0f, 0.0f, 1.0f), 0.5f * M3D_PI); // axis along X
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_CAPSULE;
	sd.radius = 0.3f;
	sd.halfHeight = 0.5f;
	m3d_body_id cap = m3d_body_create(world, &bd, &sd);

	for (int i = 0; i < 240; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	m3d_vec3 p = m3d_body_position(world, cap);
	CHECK(fabsf(p.y - 0.3f) < 0.02f, "capsule rests at y=radius");
	m3d_world_destroy(world);
}

// A 10-box tower must survive 5 simulated seconds without collapsing.
static void TestBoxStack(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	const int n = 10;
	m3d_body_id boxes[n];
	for (int i = 0; i < n; ++i)
	{
		boxes[i] = AddBox(world, m3d_v3(0.0f, 0.5f + 1.0f * (float)i, 0.0f), m3d_v3(0.5f, 0.5f, 0.5f), 0.6f);
	}

	for (int i = 0; i < 300; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}

	bool standing = true;
	for (int i = 0; i < n; ++i)
	{
		m3d_vec3 p = m3d_body_position(world, boxes[i]);
		float expectedY = 0.5f + 1.0f * (float)i;
		if (fabsf(p.y - expectedY) > 0.1f || fabsf(p.x) > 0.2f || fabsf(p.z) > 0.2f)
		{
			standing = false;
			printf("  box %d at (%.3f, %.3f, %.3f), expected y=%.1f\n", i, p.x, p.y, p.z, expectedY);
		}
	}
	CHECK(standing, "10-box stack remains standing after 5s");
	m3d_world_destroy(world);
}

// A bouncy ball must rebound to a meaningful fraction of its drop height.
static void TestRestitution(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = m3d_v3(0.0f, 2.0f, 0.0f);
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_SPHERE;
	sd.radius = 0.5f;
	sd.restitution = 0.8f;
	m3d_body_id ball = m3d_body_create(world, &bd, &sd);

	float maxHeightAfterBounce = 0.0f;
	bool bounced = false;
	for (int i = 0; i < 240; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
		m3d_vec3 p = m3d_body_position(world, ball);
		m3d_vec3 v = m3d_body_linear_velocity(world, ball);
		if (v.y > 0.5f)
		{
			bounced = true;
		}
		if (bounced && p.y > maxHeightAfterBounce)
		{
			maxHeightAfterBounce = p.y;
		}
	}
	CHECK(bounced, "restitution 0.8 ball bounces");
	// drop from 2.0 (center) => rebound center height should be well above 1.0
	CHECK(maxHeightAfterBounce > 1.0f, "bounce reaches a reasonable height");
	m3d_world_destroy(world);
}

// Bodies at rest must fall asleep; an impulse must wake them.
static void TestSleepAndWake(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);
	m3d_body_id box = AddBox(world, m3d_v3(0.0f, 0.5f, 0.0f), m3d_v3(0.5f, 0.5f, 0.5f), 0.6f);

	for (int i = 0; i < 180; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	CHECK(!m3d_body_is_awake(world, box), "resting box falls asleep");

	m3d_body_apply_impulse(world, box, m3d_v3(0.0f, 500.0f, 0.0f), m3d_body_position(world, box));
	CHECK(m3d_body_is_awake(world, box), "impulse wakes the box");
	m3d_world_step(world, 1.0f / 60.0f, 4);
	m3d_vec3 v = m3d_body_linear_velocity(world, box);
	CHECK(v.y > 0.0f, "woken box moves upward");
	m3d_world_destroy(world);
}

// Ray casting must return the closest body with a sensible normal.
static void TestRayCast(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);
	m3d_body_id box = AddBox(world, m3d_v3(0.0f, 5.0f, 0.0f), m3d_v3(0.5f, 0.5f, 0.5f), 0.6f);

	m3d_ray_result r = m3d_world_ray_cast(world, m3d_v3(0.0f, 10.0f, 0.0f), m3d_v3(0.0f, -20.0f, 0.0f));
	CHECK(r.hit, "downward ray hits");
	CHECK(r.body.index == box.index, "ray hits the box first");
	CHECK(fabsf(r.point.y - 5.5f) < 0.01f, "ray hit point on top face");
	CHECK(r.normal.y > 0.99f, "ray hit normal points up");

	m3d_ray_result miss = m3d_world_ray_cast(world, m3d_v3(100.0f, 10.0f, 100.0f), m3d_v3(0.0f, 1.0f, 0.0f));
	CHECK(!miss.hit, "ray away from all bodies misses");
	m3d_world_destroy(world);
}

// A distance joint must hold the pendulum at its length.
static void TestDistanceJoint(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);

	m3d_body_def ad = m3d_body_def_default();
	ad.type = M3D_BODY_STATIC;
	ad.position = m3d_v3(0.0f, 5.0f, 0.0f);
	m3d_shape_def asd = m3d_shape_def_default();
	asd.type = M3D_SHAPE_SPHERE;
	asd.radius = 0.1f;
	m3d_body_id anchor = m3d_body_create(world, &ad, &asd);

	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = m3d_v3(2.0f, 5.0f, 0.0f);
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_SPHERE;
	sd.radius = 0.25f;
	m3d_body_id bob = m3d_body_create(world, &bd, &sd);

	m3d_distance_joint_def jd;
	memset(&jd, 0, sizeof(jd));
	jd.bodyA = anchor;
	jd.bodyB = bob;
	jd.localAnchorA = m3d_v3_zero();
	jd.localAnchorB = m3d_v3_zero();
	jd.length = 2.0f;
	m3d_joint_create_distance(world, &jd);

	float maxErr = 0.0f;
	for (int i = 0; i < 300; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
		m3d_vec3 p = m3d_body_position(world, bob);
		float len = m3d_length(m3d_sub(p, m3d_v3(0.0f, 5.0f, 0.0f)));
		float err = fabsf(len - 2.0f);
		if (err > maxErr)
		{
			maxErr = err;
		}
	}
	CHECK(maxErr < 0.05f, "distance joint keeps pendulum length within 5cm");
	m3d_world_destroy(world);
}

// A ball joint must pin the anchor points together.
static void TestBallJoint(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);

	m3d_body_def ad = m3d_body_def_default();
	ad.type = M3D_BODY_STATIC;
	ad.position = m3d_v3(0.0f, 5.0f, 0.0f);
	m3d_shape_def asd = m3d_shape_def_default();
	asd.type = M3D_SHAPE_SPHERE;
	asd.radius = 0.1f;
	m3d_body_id anchor = m3d_body_create(world, &ad, &asd);

	// start with the pinned end exactly at the anchor
	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = m3d_v3(0.5f, 5.0f, 0.0f);
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_BOX;
	sd.halfExtents = m3d_v3(0.5f, 0.1f, 0.1f);
	m3d_body_id rod = m3d_body_create(world, &bd, &sd);

	m3d_ball_joint_def jd;
	memset(&jd, 0, sizeof(jd));
	jd.bodyA = anchor;
	jd.bodyB = rod;
	jd.localAnchorA = m3d_v3_zero();
	jd.localAnchorB = m3d_v3(-0.5f, 0.0f, 0.0f);
	m3d_joint_create_ball(world, &jd);

	float maxErr = 0.0f;
	for (int i = 0; i < 300; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
		m3d_transform xf = m3d_body_transform(world, rod);
		m3d_vec3 pinned = m3d_transform_point(xf, m3d_v3(-0.5f, 0.0f, 0.0f));
		float err = m3d_length(m3d_sub(pinned, m3d_v3(0.0f, 5.0f, 0.0f)));
		if (err > maxErr)
		{
			maxErr = err;
		}
	}
	CHECK(maxErr < 0.05f, "ball joint keeps anchors pinned within 5cm");
	m3d_world_destroy(world);
}

// Contact begin/end events fire when a body lands and when it is destroyed.
static void TestContactEvents(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);
	m3d_body_id box = AddBox(world, m3d_v3(0.0f, 2.0f, 0.0f), m3d_v3(0.5f, 0.5f, 0.5f), 0.6f);

	bool sawBegin = false;
	for (int i = 0; i < 120; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
		m3d_contact_events ev = m3d_world_contact_events(world);
		if (ev.beginCount > 0)
		{
			sawBegin = true;
		}
	}
	CHECK(sawBegin, "begin-touch event fires on landing");

	m3d_body_destroy(world, box);
	CHECK(!m3d_body_is_valid(world, box), "destroyed body id is invalid");
	m3d_world_step(world, 1.0f / 60.0f, 4);
	m3d_world_destroy(world);
}

// A fast small sphere must not tunnel through a thin static wall
// (speculative contacts).
static void TestSpeculativeNoTunnel(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);

	m3d_body_def wallDef = m3d_body_def_default();
	wallDef.type = M3D_BODY_STATIC;
	wallDef.position = m3d_v3(10.0f, 0.0f, 0.0f);
	m3d_shape_def wallShape = m3d_shape_def_default();
	wallShape.type = M3D_SHAPE_BOX;
	wallShape.halfExtents = m3d_v3(0.05f, 5.0f, 5.0f);
	m3d_body_create(world, &wallDef, &wallShape);

	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = m3d_v3(0.0f, 0.0f, 0.0f);
	bd.linearVelocity = m3d_v3(50.0f, 0.0f, 0.0f); // 0.83 m per 60Hz step
	bd.gravityScale = 0.0f;
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_SPHERE;
	sd.radius = 0.2f;
	m3d_body_id bullet = m3d_body_create(world, &bd, &sd);

	for (int i = 0; i < 60; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	m3d_vec3 p = m3d_body_position(world, bullet);
	CHECK(p.x < 10.0f, "fast sphere stopped by thin wall (no tunneling)");
	m3d_world_destroy(world);
}

// Parallel solve must produce identical results to single-threaded solve.
static void TestDeterminismAcrossThreads(void)
{
	float positions[2][3];
	for (int pass = 0; pass < 2; ++pass)
	{
		m3d_world_def wd = m3d_world_def_default();
		wd.workerCount = pass == 0 ? 1 : 4;
		m3d_world* world = m3d_world_create(&wd);
		AddGround(world);
		m3d_body_id top = { -1, 0 };
		for (int i = 0; i < 5; ++i)
		{
			for (int j = 0; j < 5; ++j)
			{
				m3d_body_id id = AddBox(world, m3d_v3(2.0f * (float)i - 4.0f, 0.5f + (float)j, 0.0f),
									  m3d_v3(0.45f, 0.45f, 0.45f), 0.6f);
				top = id;
			}
		}
		for (int i = 0; i < 120; ++i)
		{
			m3d_world_step(world, 1.0f / 60.0f, 4);
		}
		m3d_vec3 p = m3d_body_position(world, top);
		positions[pass][0] = p.x;
		positions[pass][1] = p.y;
		positions[pass][2] = p.z;
		m3d_world_destroy(world);
	}
	CHECK(positions[0][0] == positions[1][0] && positions[0][1] == positions[1][1] && positions[0][2] == positions[1][2],
		  "1-thread and 4-thread runs are bit-identical");
}

// Bodies in the same negative group must pass through each other; bodies
// with disjoint category/mask must not collide either.
static void TestCollisionFilter(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	// negative group: two overlapping-path spheres fall through each other
	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = m3d_v3(0.0f, 4.0f, 0.0f);
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_SPHERE;
	sd.radius = 0.5f;
	sd.filter.groupIndex = -7;
	m3d_body_id upper = m3d_body_create(world, &bd, &sd);
	bd.position = m3d_v3(0.0f, 1.0f, 0.0f);
	m3d_body_id lower = m3d_body_create(world, &bd, &sd);

	// disjoint masks: a "ghost" box that only collides with category 0x2
	// (the ground is category 0x1) falls straight through the ground
	bd.position = m3d_v3(5.0f, 2.0f, 0.0f);
	m3d_shape_def gh = m3d_shape_def_default();
	gh.type = M3D_SHAPE_BOX;
	gh.halfExtents = m3d_v3(0.4f, 0.4f, 0.4f);
	gh.filter.maskBits = 0x2;
	m3d_body_id ghost = m3d_body_create(world, &bd, &gh);

	for (int i = 0; i < 180; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	// both group -7 spheres rest on the ground at the same spot: they
	// ignored each other but not the ground
	m3d_vec3 pu = m3d_body_position(world, upper);
	m3d_vec3 pl = m3d_body_position(world, lower);
	CHECK(fabsf(pu.y - 0.5f) < 0.05f && fabsf(pl.y - 0.5f) < 0.05f, "same negative group: bodies interpenetrate");
	CHECK(m3d_body_position(world, ghost).y < -5.0f, "disjoint masks: body ignores the ground");
	m3d_world_destroy(world);
}

// A sensor zone must report begin/end events without exerting forces.
static void TestSensor(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	m3d_body_def zd = m3d_body_def_default();
	zd.position = m3d_v3(0.0f, 2.0f, 0.0f); // static trigger volume above the ground
	m3d_shape_def zs = m3d_shape_def_default();
	zs.type = M3D_SHAPE_BOX;
	zs.halfExtents = m3d_v3(1.0f, 1.0f, 1.0f);
	zs.isSensor = true;
	m3d_body_id zone = m3d_body_create(world, &zd, &zs);

	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = m3d_v3(0.0f, 6.0f, 0.0f);
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_SPHERE;
	sd.radius = 0.3f;
	m3d_body_id ball = m3d_body_create(world, &bd, &sd);

	int sensorBegin = 0, sensorEnd = 0;
	for (int i = 0; i < 240; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
		m3d_contact_events ev = m3d_world_contact_events(world);
		for (int k = 0; k < ev.beginCount; ++k)
		{
			if (ev.beginEvents[k].isSensor)
			{
				++sensorBegin;
			}
		}
		for (int k = 0; k < ev.endCount; ++k)
		{
			if (ev.endEvents[k].isSensor)
			{
				++sensorEnd;
			}
		}
	}
	CHECK(sensorBegin >= 1, "sensor reported a begin event");
	CHECK(sensorEnd >= 1, "sensor reported an end event after the body passed");
	m3d_vec3 p = m3d_body_position(world, ball);
	CHECK(fabsf(p.y - 0.3f) < 0.05f, "sensor exerted no force: ball fell through the zone to the ground");
	(void)zone;
	m3d_world_destroy(world);
}

// F = m*a: a constant force must produce the expected velocity, a torque
// must spin the body.
static void TestForcesAndTorques(void)
{
	m3d_world_def wd = m3d_world_def_default();
	wd.gravity = m3d_v3_zero();
	m3d_world* world = m3d_world_create(&wd);

	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.angularDamping = 0.0f;
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_BOX;
	sd.halfExtents = m3d_v3(0.5f, 0.5f, 0.5f);
	sd.density = 1000.0f; // 1 m^3 -> 1000 kg
	m3d_body_id box = m3d_body_create(world, &bd, &sd);

	for (int i = 0; i < 60; ++i)
	{
		m3d_body_apply_force(world, box, m3d_v3(2000.0f, 0.0f, 0.0f)); // a = 2 m/s^2
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	m3d_vec3 v = m3d_body_linear_velocity(world, box);
	CHECK(fabsf(v.x - 2.0f) < 0.05f, "apply_force: v = a*t after 1 s");

	for (int i = 0; i < 60; ++i)
	{
		m3d_body_apply_torque(world, box, m3d_v3(0.0f, 500.0f, 0.0f));
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	CHECK(m3d_body_angular_velocity(world, box).y > 0.5f, "apply_torque spins the body");

	// forces are cleared after each step: velocity stays constant now
	m3d_vec3 v0 = m3d_body_linear_velocity(world, box);
	m3d_world_step(world, 1.0f / 60.0f, 4);
	m3d_vec3 v1 = m3d_body_linear_velocity(world, box);
	CHECK(fabsf(v1.x - v0.x) < 1e-6f, "force does not persist past the step");
	m3d_world_destroy(world);
}

// A hinge pendulum must swing in the hinge plane only; the limit must stop
// it; the motor must spin a free wheel to the target speed.
static void TestHingeJoint(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);

	// pendulum: box hanging from a static pivot, hinge around Z
	m3d_body_def ad = m3d_body_def_default();
	ad.position = m3d_v3(0.0f, 5.0f, 0.0f);
	m3d_shape_def as = m3d_shape_def_default();
	as.type = M3D_SHAPE_BOX;
	as.halfExtents = m3d_v3(0.1f, 0.1f, 0.1f);
	m3d_body_id pivot = m3d_body_create(world, &ad, &as);

	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.position = m3d_v3(1.0f, 5.0f, 0.0f); // horizontal start: will swing down
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_BOX;
	sd.halfExtents = m3d_v3(0.4f, 0.1f, 0.1f);
	m3d_body_id arm = m3d_body_create(world, &bd, &sd);

	m3d_hinge_joint_def hd = m3d_hinge_joint_def_default();
	hd.bodyA = pivot;
	hd.bodyB = arm;
	hd.localAnchorA = m3d_v3_zero();
	hd.localAnchorB = m3d_v3(-1.0f, 0.0f, 0.0f);
	hd.localAxisA = m3d_v3(0.0f, 0.0f, 1.0f);
	m3d_joint_id hinge = m3d_joint_create_hinge(world, &hd);

	float maxZ = 0.0f;
	for (int i = 0; i < 240; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
		m3d_vec3 p = m3d_body_position(world, arm);
		float z = fabsf(p.z);
		maxZ = z > maxZ ? z : maxZ;
	}
	CHECK(maxZ < 0.02f, "hinge pendulum stays in the hinge plane");
	m3d_vec3 p = m3d_body_position(world, arm);
	float dist = m3d_length(m3d_sub(p, m3d_v3(0.0f, 5.0f, 0.0f)));
	CHECK(fabsf(dist - 1.0f) < 0.05f, "hinge point constraint holds the arm at the pivot");
	m3d_joint_destroy(world, hinge);
	m3d_body_destroy(world, arm);

	// limited hinge: pendulum released horizontally may not swing below the
	// -45 degree limit
	bd.position = m3d_v3(1.0f, 5.0f, 0.0f);
	arm = m3d_body_create(world, &bd, &sd);
	hd.bodyB = arm;
	hd.enableLimit = true;
	hd.lowerAngle = -0.25f * M3D_PI;
	hd.upperAngle = 0.25f * M3D_PI;
	hinge = m3d_joint_create_hinge(world, &hd);
	float minAngle = 0.0f;
	for (int i = 0; i < 240; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
		float a = m3d_joint_hinge_angle(world, hinge);
		minAngle = a < minAngle ? a : minAngle;
	}
	CHECK(minAngle > -0.25f * M3D_PI - 0.05f, "hinge limit holds the swing");
	m3d_joint_destroy(world, hinge);
	m3d_body_destroy(world, arm);

	// motor: a balanced wheel spun up to 5 rad/s
	bd.position = m3d_v3(0.0f, 5.0f, 0.0f);
	m3d_shape_def ws = m3d_shape_def_default();
	ws.type = M3D_SHAPE_SPHERE;
	ws.radius = 0.5f;
	m3d_body_id wheel = m3d_body_create(world, &bd, &ws);
	hd.bodyB = wheel;
	hd.localAnchorB = m3d_v3_zero();
	hd.enableLimit = false;
	hd.enableMotor = true;
	hd.motorSpeed = 5.0f;
	hd.maxMotorTorque = 1000.0f;
	hinge = m3d_joint_create_hinge(world, &hd);
	for (int i = 0; i < 120; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	CHECK(fabsf(m3d_body_angular_velocity(world, wheel).z - 5.0f) < 0.2f, "hinge motor reaches target speed");
	m3d_world_destroy(world);
}

// Welded boxes must move as one rigid body.
static void TestWeldJoint(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	m3d_body_id a = AddBox(world, m3d_v3(0.0f, 3.0f, 0.0f), m3d_v3(0.4f, 0.4f, 0.4f), 0.6f);
	m3d_body_id b = AddBox(world, m3d_v3(0.85f, 3.0f, 0.0f), m3d_v3(0.4f, 0.4f, 0.4f), 0.6f);
	m3d_weld_joint_def jd = m3d_weld_joint_def_default();
	jd.bodyA = a;
	jd.bodyB = b;
	jd.localAnchorA = m3d_v3(0.425f, 0.0f, 0.0f);
	jd.localAnchorB = m3d_v3(-0.425f, 0.0f, 0.0f);
	m3d_joint_create_weld(world, &jd);

	// drop onto the ground and keep simulating: the relative pose must hold
	for (int i = 0; i < 300; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	m3d_vec3 pa = m3d_body_position(world, a);
	m3d_vec3 pb = m3d_body_position(world, b);
	CHECK(fabsf(m3d_length(m3d_sub(pb, pa)) - 0.85f) < 0.03f, "weld keeps the anchor distance");
	m3d_quat qa = m3d_body_rotation(world, a);
	m3d_quat qb = m3d_body_rotation(world, b);
	m3d_quat rel = m3d_quat_mul(m3d_quat_conj(qa), qb);
	CHECK(fabsf(rel.w) > 0.999f, "weld keeps the relative orientation");
	m3d_world_destroy(world);
}

// Sphere cast: exact fraction against a box face, plus a clean miss.
static void TestSphereCastAndOverlap(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world); // top face at y = 0

	// downward probe from y=2 with a 0.25 sphere must hit at y=0.25 -> t=0.875
	m3d_ray_result hit = m3d_world_sphere_cast(world, m3d_v3(0.0f, 2.0f, 0.0f), 0.25f, m3d_v3(0.0f, -2.0f, 0.0f),
											   0xFFFFFFFFu);
	CHECK(hit.hit, "sphere cast hits the ground");
	CHECK(fabsf(hit.fraction - 0.875f) < 0.01f, "sphere cast fraction accounts for the radius");
	CHECK(hit.normal.y > 0.99f, "sphere cast normal points up");

	m3d_ray_result miss = m3d_world_sphere_cast(world, m3d_v3(0.0f, 5.0f, 0.0f), 0.25f, m3d_v3(0.0f, -2.0f, 0.0f),
												0xFFFFFFFFu);
	CHECK(!miss.hit, "short sphere cast misses");

	// overlap query: one box inside the region, one outside
	m3d_body_id inside = AddBox(world, m3d_v3(0.0f, 0.5f, 0.0f), m3d_v3(0.4f, 0.4f, 0.4f), 0.6f);
	AddBox(world, m3d_v3(20.0f, 0.5f, 0.0f), m3d_v3(0.4f, 0.4f, 0.4f), 0.6f);
	struct Ctx
	{
		int count;
		m3d_body_id inside;
		bool sawInside;
	} ctx = { 0, inside, false };
	m3d_aabb region = { m3d_v3(-2.0f, 0.1f, -2.0f), m3d_v3(2.0f, 2.0f, 2.0f) };
	m3d_world_overlap_aabb(world, region,
						   [](m3d_body_id id, void* c) {
							   Ctx* cx = (Ctx*)c;
							   ++cx->count;
							   if (id.index == cx->inside.index)
							   {
								   cx->sawInside = true;
							   }
							   return true;
						   },
						   &ctx);
	CHECK(ctx.sawInside, "overlap query finds the box in the region");
	CHECK(ctx.count == 1, "overlap query skips bodies outside the region");
	m3d_world_destroy(world);
}

// A very fast small body must not pass through a thin wall: per-substep
// travel exceeds the wall span, so the fresh manifold flips to the far
// face - the tunneling rescue must keep the pre-crossing manifold and
// push the body back out the side it entered from.
static void TestThinWallHighSpeed(void)
{
	m3d_world_def wd = m3d_world_def_default();
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	m3d_body_def wallD = m3d_body_def_default();
	wallD.position = m3d_v3(20.0f, 2.0f, 0.0f);
	m3d_shape_def wallS = m3d_shape_def_default();
	wallS.type = M3D_SHAPE_BOX;
	wallS.halfExtents = m3d_v3(0.025f, 2.0f, 2.0f); // 5 cm wall
	m3d_body_create(world, &wallD, &wallS);

	m3d_body_def peaD = m3d_body_def_default();
	peaD.type = M3D_BODY_DYNAMIC;
	peaD.position = m3d_v3(12.0f, 2.0f, 0.0f);
	peaD.linearVelocity = m3d_v3(150.0f, 0.0f, 0.0f); // 0.625 m per substep
	peaD.gravityScale = 0.0f;
	m3d_shape_def peaS = m3d_shape_def_default();
	peaS.type = M3D_SHAPE_SPHERE;
	peaS.radius = 0.06f;
	m3d_body_id pea = m3d_body_create(world, &peaD, &peaS);

	for (int i = 0; i < 60; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, 4);
	}
	m3d_vec3 p = m3d_body_position(world, pea);
	CHECK(p.x < 20.0f, "150 m/s pea ends on the near side of a 5 cm wall");
	m3d_world_destroy(world);
}

// Joints must not break determinism across thread counts.
static void TestJointDeterminism(void)
{
	float positions[2][3];
	for (int pass = 0; pass < 2; ++pass)
	{
		m3d_world_def wd = m3d_world_def_default();
		wd.workerCount = pass == 0 ? 1 : 4;
		m3d_world* world = m3d_world_create(&wd);
		AddGround(world);

		// chain of hinged boxes swinging from a static pivot
		m3d_body_def ad = m3d_body_def_default();
		ad.position = m3d_v3(0.0f, 8.0f, 0.0f);
		m3d_shape_def as = m3d_shape_def_default();
		as.type = M3D_SHAPE_BOX;
		as.halfExtents = m3d_v3(0.1f, 0.1f, 0.1f);
		m3d_body_id prev = m3d_body_create(world, &ad, &as);
		m3d_body_id last = prev;
		for (int i = 0; i < 6; ++i)
		{
			m3d_body_id link = AddBox(world, m3d_v3(0.9f * (float)(i + 1), 8.0f, 0.0f), m3d_v3(0.4f, 0.15f, 0.15f), 0.6f);
			m3d_hinge_joint_def hd = m3d_hinge_joint_def_default();
			hd.bodyA = prev;
			hd.bodyB = link;
			hd.localAnchorA = i == 0 ? m3d_v3_zero() : m3d_v3(0.45f, 0.0f, 0.0f);
			hd.localAnchorB = m3d_v3(-0.45f, 0.0f, 0.0f);
			hd.localAxisA = m3d_v3(0.0f, 0.0f, 1.0f);
			m3d_joint_create_hinge(world, &hd);
			prev = link;
			last = link;
		}
		// some boxes raining on the chain to couple contacts and joints
		for (int i = 0; i < 12; ++i)
		{
			AddBox(world, m3d_v3(0.6f * (float)i - 2.0f, 10.0f + 0.3f * (float)i, 0.05f * (float)i),
				   m3d_v3(0.25f, 0.25f, 0.25f), 0.6f);
		}
		for (int i = 0; i < 180; ++i)
		{
			m3d_world_step(world, 1.0f / 60.0f, 4);
		}
		m3d_vec3 p = m3d_body_position(world, last);
		positions[pass][0] = p.x;
		positions[pass][1] = p.y;
		positions[pass][2] = p.z;
		m3d_world_destroy(world);
	}
	CHECK(positions[0][0] == positions[1][0] && positions[0][1] == positions[1][1] &&
			  positions[0][2] == positions[1][2],
		  "hinge chain is bit-identical across thread counts");
}

// Real-world friction laws on an inclined plane: static hold below the
// cone, Coulomb sliding above it with a = g (sin th - mu cos th).
static void TestInclinedPlane(void)
{
	// shared incline builder
	auto scene = [](m3d_world** outW, float thetaDeg, float mu) {
		float th = thetaDeg * M3D_PI / 180.0f;
		m3d_world_def wd = m3d_world_def_default();
		m3d_world* w = m3d_world_create(&wd);
		m3d_body_def gd = m3d_body_def_default();
		gd.orientation = m3d_quat_axis_angle(m3d_v3(0, 0, 1), th);
		m3d_shape_def gs = m3d_shape_def_default();
		gs.type = M3D_SHAPE_BOX;
		gs.halfExtents = m3d_v3(60, 1, 8);
		gs.friction = mu;
		m3d_body_create(w, &gd, &gs);
		m3d_vec3 nrm = m3d_v3(-sinf(th), cosf(th), 0);
		m3d_vec3 along = m3d_v3(cosf(th), sinf(th), 0);
		m3d_body_def bd = m3d_body_def_default();
		bd.type = M3D_BODY_DYNAMIC;
		bd.position = m3d_add(m3d_scale(20.0f, along), m3d_scale(1.402f, nrm));
		bd.orientation = gd.orientation;
		m3d_shape_def sd = m3d_shape_def_default();
		sd.type = M3D_SHAPE_BOX;
		sd.halfExtents = m3d_v3(0.4f, 0.4f, 0.4f);
		sd.friction = mu;
		*outW = w;
		return m3d_body_create(w, &bd, &sd);
	};

	// 28 deg, mu 0.6 (tan 28 = 0.53): static friction holds
	m3d_world* w;
	m3d_body_id box = scene(&w, 28.0f, 0.6f);
	m3d_vec3 p0 = m3d_body_position(w, box);
	for (int i = 0; i < 480; ++i)
	{
		m3d_world_step(w, 1.0f / 60.0f, 4);
	}
	CHECK(m3d_length(m3d_sub(m3d_body_position(w, box), p0)) < 0.05f, "box holds below the friction cone");
	m3d_world_destroy(w);

	// 38 deg, mu 0.6 (tan 38 = 0.78, past even the ~1.15x static bonus): slides
	box = scene(&w, 38.0f, 0.6f);
	p0 = m3d_body_position(w, box);
	for (int i = 0; i < 480; ++i)
	{
		m3d_world_step(w, 1.0f / 60.0f, 4);
	}
	CHECK(m3d_length(m3d_sub(m3d_body_position(w, box), p0)) > 0.5f, "box slides above the friction cone");
	m3d_world_destroy(w);

	// kinetic: 35 deg, mu 0.3 -> a = g (sin - mu cos) = 3.22 m/s^2
	box = scene(&w, 35.0f, 0.3f);
	for (int i = 0; i < 60; ++i)
	{
		m3d_world_step(w, 1.0f / 240.0f, 1); // break loose
	}
	m3d_vec3 v1 = m3d_body_linear_velocity(w, box);
	for (int i = 0; i < 240; ++i)
	{
		m3d_world_step(w, 1.0f / 240.0f, 1);
	}
	float aMeas = m3d_length(m3d_sub(m3d_body_linear_velocity(w, box), v1));
	float th = 35.0f * M3D_PI / 180.0f;
	float aExp = 9.81f * (sinf(th) - 0.3f * cosf(th));
	CHECK(fabsf(aMeas - aExp) < 0.1f * aExp, "Coulomb sliding acceleration matches g(sin-mu cos)");
	m3d_world_destroy(w);
}

// Rolling dynamics: a solid sphere rolls down at 5/7 g sin(th), converts
// backspin to translation at v = 2/7 w r, and a free roller must NOT
// self-accelerate (regression for the contact-arm torque pump).
static void TestRollingPhysics(void)
{
	// rolling downhill
	{
		float th = 20.0f * M3D_PI / 180.0f;
		m3d_world_def wd = m3d_world_def_default();
		m3d_world* w = m3d_world_create(&wd);
		m3d_body_def gd = m3d_body_def_default();
		gd.orientation = m3d_quat_axis_angle(m3d_v3(0, 0, 1), th);
		m3d_shape_def gs = m3d_shape_def_default();
		gs.type = M3D_SHAPE_BOX;
		gs.halfExtents = m3d_v3(60, 1, 8);
		gs.friction = 0.8f;
		m3d_body_create(w, &gd, &gs);
		m3d_body_def bd = m3d_body_def_default();
		bd.type = M3D_BODY_DYNAMIC;
		bd.position = m3d_add(m3d_scale(20.0f, m3d_v3(cosf(th), sinf(th), 0)),
							  m3d_scale(1.402f, m3d_v3(-sinf(th), cosf(th), 0)));
		bd.angularDamping = 0.0f;
		m3d_shape_def sd = m3d_shape_def_default();
		sd.type = M3D_SHAPE_SPHERE;
		sd.radius = 0.4f;
		sd.friction = 0.8f;
		sd.rollingResistance = 0.0f;
		m3d_body_id ball = m3d_body_create(w, &bd, &sd);
		for (int i = 0; i < 120; ++i)
		{
			m3d_world_step(w, 1.0f / 240.0f, 1);
		}
		m3d_vec3 v1 = m3d_body_linear_velocity(w, ball);
		for (int i = 0; i < 240; ++i)
		{
			m3d_world_step(w, 1.0f / 240.0f, 1);
		}
		float aMeas = m3d_length(m3d_sub(m3d_body_linear_velocity(w, ball), v1));
		float aExp = 5.0f / 7.0f * 9.81f * sinf(th);
		CHECK(fabsf(aMeas - aExp) < 0.08f * aExp, "solid sphere rolls at 5/7 g sin(th)");
	m3d_world_destroy(w);
	}
	// backspin conversion + no self-acceleration
	{
		m3d_world_def wd = m3d_world_def_default();
		m3d_world* w = m3d_world_create(&wd);
		AddGround(w);
		m3d_body_def bd = m3d_body_def_default();
		bd.type = M3D_BODY_DYNAMIC;
		bd.position = m3d_v3(0, 0.5f, 0);
		bd.angularVelocity = m3d_v3(0, 0, -10.0f);
		bd.angularDamping = 0.0f;
		m3d_shape_def sd = m3d_shape_def_default();
		sd.type = M3D_SHAPE_SPHERE;
		sd.radius = 0.5f;
		sd.friction = 0.5f;
		sd.rollingResistance = 0.0f;
		m3d_body_id ball = m3d_body_create(w, &bd, &sd);
		for (int i = 0; i < 150; ++i)
		{
			m3d_world_step(w, 1.0f / 240.0f, 1); // spin-up completes (~0.3 s)
		}
		float vRoll = m3d_body_linear_velocity(w, ball).x;
		CHECK(fabsf(vRoll - 2.0f / 7.0f * 10.0f * 0.5f) < 0.1f, "backspin converts at v = 2/7 w r");
		for (int i = 0; i < 480; ++i)
		{
			m3d_world_step(w, 1.0f / 240.0f, 1);
		}
		float vLater = m3d_body_linear_velocity(w, ball).x;
		CHECK(vLater <= vRoll + 0.02f, "free roller never self-accelerates");
		m3d_world_destroy(w);
	}
}

// The intermediate-axis (Dzhanibekov) instability requires the gyroscopic
// term: a brick spun about its middle axis must flip.
static void TestDzhanibekov(void)
{
	m3d_world_def wd = m3d_world_def_default();
	wd.gravity = m3d_v3_zero();
	m3d_world* w = m3d_world_create(&wd);
	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_DYNAMIC;
	bd.angularDamping = 0.0f;
	bd.angularVelocity = m3d_v3(0.02f, 8.0f, 0.0f);
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_BOX;
	sd.halfExtents = m3d_v3(0.5f, 0.3f, 0.1f);
	m3d_body_id brick = m3d_body_create(w, &bd, &sd);
	float minDot = 1.0f;
	for (int i = 0; i < 2400; ++i)
	{
		m3d_world_step(w, 1.0f / 240.0f, 1);
		float d = m3d_rotate(m3d_body_rotation(w, brick), m3d_v3(0, 1, 0)).y;
		minDot = d < minDot ? d : minDot;
	}
	CHECK(minDot < -0.5f, "Dzhanibekov flip about the intermediate axis");
	m3d_world_destroy(w);
}

int main(void)
{
	TestSphereRest();
	TestCapsuleRest();
	TestBoxStack();
	TestRestitution();
	TestSleepAndWake();
	TestRayCast();
	TestDistanceJoint();
	TestBallJoint();
	TestContactEvents();
	TestSpeculativeNoTunnel();
	TestDeterminismAcrossThreads();
	TestCollisionFilter();
	TestSensor();
	TestForcesAndTorques();
	TestHingeJoint();
	TestWeldJoint();
	TestSphereCastAndOverlap();
	TestThinWallHighSpeed();
	TestInclinedPlane();
	TestRollingPhysics();
	TestDzhanibekov();
	TestJointDeterminism();

	printf("\n%d checks, %d failures\n", g_testCount, g_failCount);
	return g_failCount == 0 ? 0 : 1;
}
