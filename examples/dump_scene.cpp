// Runs a showcase scene and dumps per-frame body transforms to a binary
// file for the Python renderer. Bodies are all created up front so the
// frame layout is fixed.
//
// The scene: a wrecking ball on a chain of ball joints swings down and
// demolishes a box pyramid, then a rain of mixed shapes settles over the
// debris. Joints, stacking, impacts and mixed-shape piles in one shot.
#include "melon3d.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>

struct BodyInfo
{
	m3d_body_id id;
	int32_t shapeType;
	float r, hh, ex, ey, ez;
	// 0 ground, 1 pyramid, 2 wrecking ball, 3 rain-sphere, 4 rain-capsule,
	// 5 rain-box, 6 chain link, 7 anchor pillar
	int32_t colorId;
	int32_t isStatic;
};

int main()
{
	m3d_world_def wd = m3d_world_def_default();
	wd.workerCount = 8;
	m3d_world* world = m3d_world_create(&wd);

	std::vector<BodyInfo> infos;

	auto addBox = [&](m3d_vec3 p, m3d_vec3 he, int32_t color, bool isStatic, m3d_quat q) {
		m3d_body_def bd = m3d_body_def_default();
		bd.type = isStatic ? M3D_BODY_STATIC : M3D_BODY_DYNAMIC;
		bd.position = p;
		bd.orientation = q;
		m3d_shape_def sd = m3d_shape_def_default();
		sd.type = M3D_SHAPE_BOX;
		sd.halfExtents = he;
		sd.friction = 0.7f;
		BodyInfo bi;
		bi.id = m3d_body_create(world, &bd, &sd);
		bi.shapeType = M3D_SHAPE_BOX;
		bi.r = 0; bi.hh = 0; bi.ex = he.x; bi.ey = he.y; bi.ez = he.z;
		bi.colorId = color; bi.isStatic = isStatic ? 1 : 0;
		infos.push_back(bi);
		return bi.id;
	};
	auto addSphere = [&](m3d_vec3 p, float r, int32_t color, float density, m3d_vec3 v) {
		m3d_body_def bd = m3d_body_def_default();
		bd.type = M3D_BODY_DYNAMIC;
		bd.position = p;
		bd.linearVelocity = v;
		m3d_shape_def sd = m3d_shape_def_default();
		sd.type = M3D_SHAPE_SPHERE;
		sd.radius = r;
		sd.friction = 0.5f;
		sd.restitution = 0.15f;
		sd.density = density;
		BodyInfo bi;
		bi.id = m3d_body_create(world, &bd, &sd);
		bi.shapeType = M3D_SHAPE_SPHERE;
		bi.r = r; bi.hh = 0; bi.ex = bi.ey = bi.ez = 0;
		bi.colorId = color; bi.isStatic = 0;
		infos.push_back(bi);
		return bi.id;
	};
	auto addCapsule = [&](m3d_vec3 p, float r, float hh, int32_t color, m3d_quat q, float density, m3d_vec3 v) {
		m3d_body_def bd = m3d_body_def_default();
		bd.type = M3D_BODY_DYNAMIC;
		bd.position = p;
		bd.orientation = q;
		bd.linearVelocity = v;
		m3d_shape_def sd = m3d_shape_def_default();
		sd.type = M3D_SHAPE_CAPSULE;
		sd.radius = r;
		sd.halfHeight = hh;
		sd.friction = 0.6f;
		sd.density = density;
		BodyInfo bi;
		bi.id = m3d_body_create(world, &bd, &sd);
		bi.shapeType = M3D_SHAPE_CAPSULE;
		bi.r = r; bi.hh = hh; bi.ex = bi.ey = bi.ez = 0;
		bi.colorId = color; bi.isStatic = 0;
		infos.push_back(bi);
		return bi.id;
	};

	// ground
	addBox(m3d_v3(0, -1, 0), m3d_v3(60, 1, 60), 0, true, m3d_quat_identity());

	// pyramid of boxes (base 6) at the origin
	const float he = 0.5f;
	const int base = 6;
	for (int row = 0; row < base; ++row)
	{
		int n = base - row;
		float y = he + 2.0f * he * row;
		float x0 = -he * (n - 1);
		for (int i = 0; i < n; ++i)
		{
			m3d_body_def bd = m3d_body_def_default();
			bd.type = M3D_BODY_DYNAMIC;
			bd.position = m3d_v3(x0 + 2.0f * he * i, y, 0);
			m3d_shape_def sd = m3d_shape_def_default();
			sd.type = M3D_SHAPE_BOX;
			sd.halfExtents = m3d_v3(he, he, he);
			sd.friction = 0.72f;
			BodyInfo bi;
			bi.id = m3d_body_create(world, &bd, &sd);
			bi.shapeType = M3D_SHAPE_BOX;
			bi.r = 0; bi.hh = 0; bi.ex = he; bi.ey = he; bi.ez = he;
			bi.colorId = 1; bi.isStatic = 0;
			infos.push_back(bi);
		}
	}

	// wrecking rig: a static gantry (two pillars + crossbeam, offset in Z so
	// the swing plane stays clear), a chain of capsule links on ball joints,
	// and a heavy ball. A rigid distance joint from the beam to the ball acts
	// as the pendulum rod - it pins the ball to the exact swing arc while the
	// chain provides the visible slack dynamics.
	{
		// gantry sits at x=-2 so the bottom of the swing arc is PAST the
		// pyramid center: the ball plows through with follow-through instead
		// of stopping at the middle
		const m3d_vec3 anchorPos = m3d_v3(-2.0f, 11.1f, 0.0f);
		const int linkCount = 6;
		const float linkR = 0.18f, linkHH = 0.55f;
		const float linkLen = 2.0f * (linkHH + linkR); // joint-to-joint
		const float ballR = 1.1f;
		const float theta = 1.25f; // radians from vertical at release

		addBox(m3d_v3(-2.0f, 5.6f, -2.6f), m3d_v3(0.28f, 5.6f, 0.28f), 7, true, m3d_quat_identity());
		addBox(m3d_v3(-2.0f, 5.6f, 2.6f), m3d_v3(0.28f, 5.6f, 0.28f), 7, true, m3d_quat_identity());
		m3d_body_id beam = addBox(m3d_v3(-2.0f, 11.4f, 0), m3d_v3(0.3f, 0.3f, 2.9f), 7, true, m3d_quat_identity());

		// chain direction at release (in the swing plane, +X up-tilt)
		m3d_vec3 dir = m3d_v3(sinf(theta), -cosf(theta), 0.0f);
		// capsule local +Y must map onto -dir (top of the link points at the
		// anchor): Rz(theta) * (0,1,0) = (-sin, cos, 0) = -dir
		m3d_quat q = m3d_quat_axis_angle(m3d_v3(0, 0, 1), theta);

		// released with some angular speed already (a crane would not drop it
		// from rest): every pendulum body starts with v = omega x r
		const m3d_vec3 omega = m3d_v3(0, 0, -0.55f);
		auto swingVel = [&](m3d_vec3 c) { return m3d_cross(omega, m3d_sub(c, anchorPos)); };

		m3d_body_id prev = beam;
		for (int i = 0; i < linkCount; ++i)
		{
			m3d_vec3 c = m3d_add(anchorPos, m3d_scale((0.5f + (float)i) * linkLen, dir));
			m3d_body_id link = addCapsule(c, linkR, linkHH, 6, q, 7800.0f, swingVel(c));
			m3d_ball_joint_def jd;
			jd.bodyA = prev;
			jd.bodyB = link;
			jd.localAnchorA = i == 0 ? m3d_v3(0, -0.3f, 0) : m3d_v3(0, -(linkHH + linkR), 0);
			jd.localAnchorB = m3d_v3(0, linkHH + linkR, 0);
			m3d_joint_create_ball(world, &jd);
			prev = link;
		}
		float chainLen = linkCount * linkLen;
		m3d_vec3 ballC = m3d_add(anchorPos, m3d_scale(chainLen + ballR * 0.85f, dir));
		// heavy: the rigid rod carries the arc, so the chain does not fold
		m3d_body_id ball = addSphere(ballC, ballR, 2, 5500.0f, swingVel(ballC));
		m3d_ball_joint_def jd;
		jd.bodyA = prev;
		jd.bodyB = ball;
		jd.localAnchorA = m3d_v3(0, -(linkHH + linkR), 0);
		jd.localAnchorB = m3d_v3(0, ballR * 0.85f, 0);
		m3d_joint_create_ball(world, &jd);

		// the invisible pendulum rod: beam underside to ball center, rigid
		m3d_distance_joint_def rod;
		rod.bodyA = beam;
		rod.bodyB = ball;
		rod.localAnchorA = m3d_v3(0, -0.3f, 0);
		rod.localAnchorB = m3d_v3_zero();
		rod.length = 0.0f; // rigid at the creation distance
		rod.hertz = 0.0f;
		rod.dampingRatio = 0.0f;
		m3d_joint_create_distance(world, &rod);
	}

	// rain of mixed shapes, staggered high: the demolition reads first, then
	// they fall in and settle over the debris
	unsigned seed = 20260708u;
	auto frand = [&]() { seed = seed * 1664525u + 1013904223u; return (float)(seed >> 8) / (float)0x00FFFFFF; };
	for (int i = 0; i < 46; ++i)
	{
		float x = 9.0f * frand() - 4.5f;
		float z = 6.0f * frand() - 3.0f;
		float y = 19.0f + 2.0f * i; // higher + more stagger => later rain
		int kind = i % 3;
		if (kind == 0)
			addSphere(m3d_v3(x, y, z), 0.32f + 0.15f * frand(), 3, 1000.0f, m3d_v3_zero());
		else if (kind == 1)
			addCapsule(m3d_v3(x, y, z), 0.24f, 0.34f, 4,
					   m3d_quat_axis_angle(m3d_v3(frand(), frand(), frand()), frand() * 6.28f), 1000.0f,
					   m3d_v3_zero());
		else
			addBox(m3d_v3(x, y, z), m3d_v3(0.33f, 0.33f, 0.33f), 5, false,
				   m3d_quat_axis_angle(m3d_v3(frand(), frand(), frand()), frand() * 6.28f));
	}

	const int32_t numBodies = (int32_t)infos.size();
	const int32_t numFrames = 500;

	FILE* f = fopen("scene.bin", "wb");
	fwrite(&numFrames, sizeof(int32_t), 1, f);
	fwrite(&numBodies, sizeof(int32_t), 1, f);
	for (const BodyInfo& bi : infos)
	{
		fwrite(&bi.shapeType, sizeof(int32_t), 1, f);
		fwrite(&bi.r, sizeof(float), 1, f);
		fwrite(&bi.hh, sizeof(float), 1, f);
		fwrite(&bi.ex, sizeof(float), 1, f);
		fwrite(&bi.ey, sizeof(float), 1, f);
		fwrite(&bi.ez, sizeof(float), 1, f);
		fwrite(&bi.colorId, sizeof(int32_t), 1, f);
		fwrite(&bi.isStatic, sizeof(int32_t), 1, f);
	}

	for (int frame = 0; frame < numFrames; ++frame)
	{
		// 8 substeps: offline render, so buy extra chain stiffness for free
		m3d_world_step(world, 1.0f / 60.0f, 8);
		for (const BodyInfo& bi : infos)
		{
			m3d_transform xf = m3d_body_transform(world, bi.id);
			float d[7] = { xf.p.x, xf.p.y, xf.p.z, xf.q.x, xf.q.y, xf.q.z, xf.q.w };
			fwrite(d, sizeof(float), 7, f);
		}
	}
	fclose(f);

	m3d_world_destroy(world);
	printf("wrote scene.bin: %d frames, %d bodies\n", numFrames, numBodies);
	return 0;
}
