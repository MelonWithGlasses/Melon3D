// Runs a showcase scene and dumps per-frame body transforms to a binary
// file for the Python renderer. Bodies are all created up front so the
// frame layout is fixed.
#include "melon3d.h"
#include <stdint.h>
#include <stdio.h>
#include <vector>

struct BodyInfo
{
	m3d_body_id id;
	int32_t shapeType;
	float r, hh, ex, ey, ez;
	int32_t colorId; // 0 ground, 1 pyramid, 2 projectile, 3 rain-sphere, 4 rain-capsule, 5 rain-box
	int32_t isStatic;
};

int main()
{
	m3d_world_def wd = m3d_world_def_default();
	wd.workerCount = 8;
	m3d_world* world = m3d_world_create(&wd);

	std::vector<BodyInfo> infos;

	auto addBox = [&](m3d_vec3 p, m3d_vec3 he, int32_t color, bool isStatic, m3d_vec3 v) {
		m3d_body_def bd = m3d_body_def_default();
		bd.type = isStatic ? M3D_BODY_STATIC : M3D_BODY_DYNAMIC;
		bd.position = p;
		bd.linearVelocity = v;
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
	};
	auto addSphere = [&](m3d_vec3 p, float r, int32_t color, m3d_vec3 v, float gscale) {
		m3d_body_def bd = m3d_body_def_default();
		bd.type = M3D_BODY_DYNAMIC;
		bd.position = p;
		bd.linearVelocity = v;
		bd.gravityScale = gscale;
		m3d_shape_def sd = m3d_shape_def_default();
		sd.type = M3D_SHAPE_SPHERE;
		sd.radius = r;
		sd.friction = 0.5f;
		sd.restitution = 0.2f;
		BodyInfo bi;
		bi.id = m3d_body_create(world, &bd, &sd);
		bi.shapeType = M3D_SHAPE_SPHERE;
		bi.r = r; bi.hh = 0; bi.ex = bi.ey = bi.ez = 0;
		bi.colorId = color; bi.isStatic = 0;
		infos.push_back(bi);
	};
	auto addCapsule = [&](m3d_vec3 p, float r, float hh, int32_t color, m3d_quat q) {
		m3d_body_def bd = m3d_body_def_default();
		bd.type = M3D_BODY_DYNAMIC;
		bd.position = p;
		bd.orientation = q;
		m3d_shape_def sd = m3d_shape_def_default();
		sd.type = M3D_SHAPE_CAPSULE;
		sd.radius = r;
		sd.halfHeight = hh;
		sd.friction = 0.6f;
		BodyInfo bi;
		bi.id = m3d_body_create(world, &bd, &sd);
		bi.shapeType = M3D_SHAPE_CAPSULE;
		bi.r = r; bi.hh = hh; bi.ex = bi.ey = bi.ez = 0;
		bi.colorId = color; bi.isStatic = 0;
		infos.push_back(bi);
	};

	// ground
	addBox(m3d_v3(0, -1, 0), m3d_v3(60, 1, 60), 0, true, m3d_v3_zero());

	// pyramid of boxes (base 6)
	const float he = 0.5f;
	const int base = 6;
	for (int row = 0; row < base; ++row)
	{
		int n = base - row;
		float y = he + 2.0f * he * row;
		float x0 = -he * (n - 1);
		for (int i = 0; i < n; ++i)
		{
			addBox(m3d_v3(x0 + 2.0f * he * i, y, 0), m3d_v3(he, he, he), 1, false, m3d_v3_zero());
		}
	}

	// projectile: heavy fast sphere aimed at the mid pyramid to punch through
	addSphere(m3d_v3(-34, 3.2f, 0), 0.9f, 2, m3d_v3(23.0f, 0.8f, 0), 0.12f);

	// rain of mixed shapes, staggered high so the smash reads first, then
	// they fall in and settle over the debris
	unsigned seed = 20260705u;
	auto frand = [&]() { seed = seed * 1664525u + 1013904223u; return (float)(seed >> 8) / (float)0x00FFFFFF; };
	for (int i = 0; i < 48; ++i)
	{
		float x = 9.0f * frand() - 4.5f;
		float z = 6.0f * frand() - 3.0f;
		float y = 17.0f + 1.9f * i; // higher + more stagger => later rain
		int kind = i % 3;
		if (kind == 0)
			addSphere(m3d_v3(x, y, z), 0.35f + 0.15f * frand(), 3, m3d_v3_zero(), 1.0f);
		else if (kind == 1)
			addCapsule(m3d_v3(x, y, z), 0.25f, 0.35f, 4, m3d_quat_axis_angle(m3d_v3(frand(), frand(), frand()), frand() * 6.28f));
		else
			addBox(m3d_v3(x, y, z), m3d_v3(0.35f, 0.35f, 0.35f), 5, false, m3d_v3_zero());
	}

	const int32_t numBodies = (int32_t)infos.size();
	const int32_t numFrames = 470;

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
		m3d_world_step(world, 1.0f / 60.0f, 4);
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
