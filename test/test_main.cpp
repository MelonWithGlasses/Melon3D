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

	printf("\n%d checks, %d failures\n", g_testCount, g_failCount);
	return g_failCount == 0 ? 0 : 1;
}
