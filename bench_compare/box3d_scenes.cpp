// Benchmark scenes for erincatto/box3d, mirroring melon3d/bench/bench.cpp:
// same geometry, same 60 Hz timestep, same step counts. Timed externally
// around b3World_Step. Single-threaded (no task system callbacks provided).
#include "box3d/box3d.h"

#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using Clock = std::chrono::steady_clock;

static b3WorldId CreateWorldStd()
{
	b3WorldDef wd = b3DefaultWorldDef();
	wd.gravity = { 0.0f, -9.81f, 0.0f };
	return b3CreateWorld(&wd);
}

static void AddGround(b3WorldId world)
{
	b3BodyDef bd = b3DefaultBodyDef();
	bd.position = { 0.0, -1.0, 0.0 };
	b3BodyId body = b3CreateBody(world, &bd);
	b3ShapeDef sd = b3DefaultShapeDef();
	sd.density = 1000.0f;
	sd.baseMaterial.friction = 0.6f;
	b3BoxHull hull = b3MakeBoxHull(200.0f, 1.0f, 200.0f);
	b3CreateHullShape(body, &sd, &hull.base);
}

static void AddBox(b3WorldId world, double x, double y, double z, float he)
{
	b3BodyDef bd = b3DefaultBodyDef();
	bd.type = b3_dynamicBody;
	bd.position = { x, y, z };
	b3BodyId body = b3CreateBody(world, &bd);
	b3ShapeDef sd = b3DefaultShapeDef();
	sd.density = 1000.0f;
	sd.baseMaterial.friction = 0.6f;
	b3BoxHull hull = b3MakeCubeHull(he);
	b3CreateHullShape(body, &sd, &hull.base);
}

static void RunScene(const char* name, b3WorldId world, int steps, int substeps)
{
	double totalMs = 0.0;
	double maxMs = 0.0;
	for (int i = 0; i < steps; ++i)
	{
		auto t0 = Clock::now();
		b3World_Step(world, 1.0f / 60.0f, substeps);
		auto t1 = Clock::now();
		double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
		totalMs += ms;
		if (ms > maxMs)
		{
			maxMs = ms;
		}
	}
	printf("%-24s | avg %7.3f ms  max %7.3f ms\n", name, totalMs / steps, maxMs);
}

static void BenchPyramid(int base, int substeps)
{
	b3WorldId world = CreateWorldStd();
	AddGround(world);

	const float he = 0.5f;
	for (int row = 0; row < base; ++row)
	{
		int countInRow = base - row;
		double y = he + 2.0f * he * (double)row;
		double x0 = -he * (double)(countInRow - 1);
		for (int i = 0; i < countInRow; ++i)
		{
			AddBox(world, x0 + 2.0f * he * (double)i, y, 0.0, he);
		}
	}

	char name[64];
	snprintf(name, sizeof(name), "pyramid(%d) sub=%d", base, substeps);
	RunScene(name, world, 240, substeps);
	b3DestroyWorld(world);
}

static void BenchStackGrid(int grid, int height, int substeps)
{
	b3WorldId world = CreateWorldStd();
	AddGround(world);

	for (int gx = 0; gx < grid; ++gx)
	{
		for (int gz = 0; gz < grid; ++gz)
		{
			double x = 4.0 * gx - 2.0 * grid;
			double z = 4.0 * gz - 2.0 * grid;
			for (int j = 0; j < height; ++j)
			{
				AddBox(world, x, 0.5 + 1.0 * j, z, 0.5f);
			}
		}
	}

	char name[64];
	snprintf(name, sizeof(name), "stacks %dx%dx%d sub=%d", grid, grid, height, substeps);
	RunScene(name, world, 240, substeps);
	b3DestroyWorld(world);
}

static void BenchMixedRain(int count, int substeps)
{
	b3WorldId world = CreateWorldStd();
	AddGround(world);

	unsigned seed = 12345;
	auto frand = [&seed]() {
		seed = seed * 1664525u + 1013904223u;
		return (float)(seed >> 8) / (float)0x00FFFFFF;
	};

	for (int i = 0; i < count; ++i)
	{
		b3BodyDef bd = b3DefaultBodyDef();
		bd.type = b3_dynamicBody;
		bd.position = { 20.0f * frand() - 10.0f, 1.0 + 0.1 * i, 20.0f * frand() - 10.0f };
		b3BodyId body = b3CreateBody(world, &bd);
		b3ShapeDef sd = b3DefaultShapeDef();
		sd.density = 1000.0f;
		sd.baseMaterial.friction = 0.6f;
		int kind = i % 3;
		if (kind == 0)
		{
			b3Sphere s = { { 0.0f, 0.0f, 0.0f }, 0.3f };
			b3CreateSphereShape(body, &sd, &s);
		}
		else if (kind == 1)
		{
			b3Capsule c = { { 0.0f, -0.3f, 0.0f }, { 0.0f, 0.3f, 0.0f }, 0.2f };
			b3CreateCapsuleShape(body, &sd, &c);
		}
		else
		{
			b3BoxHull hull = b3MakeCubeHull(0.3f);
			b3CreateHullShape(body, &sd, &hull.base);
		}
	}

	char name[64];
	snprintf(name, sizeof(name), "mixed rain %d sub=%d", count, substeps);
	RunScene(name, world, 300, substeps);
	b3DestroyWorld(world);
}

static void BenchMiniTowers(int grid, int height, int substeps)
{
	b3WorldId world = CreateWorldStd();
	AddGround(world);

	for (int gx = 0; gx < grid; ++gx)
	{
		for (int gz = 0; gz < grid; ++gz)
		{
			double x = 3.0 * gx - 1.5 * grid;
			double z = 3.0 * gz - 1.5 * grid;
			for (int j = 0; j < height; ++j)
			{
				AddBox(world, x, 0.5 + 1.0 * j, z, 0.5f);
			}
		}
	}

	char name[64];
	snprintf(name, sizeof(name), "towers %dx%dx%d sub=%d", grid, grid, height, substeps);
	RunScene(name, world, 240, substeps);
	b3DestroyWorld(world);
}

// Streaming churn: bodies continuously destroyed and respawned while
// raining down. Create/destroy time is part of the measurement.
static void BenchChurn(int count, int churnPerStep, int substeps)
{
	b3WorldId world = CreateWorldStd();
	AddGround(world);

	unsigned seed = 777;
	auto frand = [&seed]() {
		seed = seed * 1664525u + 1013904223u;
		return (float)(seed >> 8) / (float)0x00FFFFFF;
	};

	std::vector<b3BodyId> ring((size_t)count);
	auto spawn = [&](int slot) {
		b3BodyDef bd = b3DefaultBodyDef();
		bd.type = b3_dynamicBody;
		bd.position = { 24.0f * frand() - 12.0f, 6.0f + 8.0f * frand(), 24.0f * frand() - 12.0f };
		b3BodyId body = b3CreateBody(world, &bd);
		b3ShapeDef sd = b3DefaultShapeDef();
		sd.density = 1000.0f;
		sd.baseMaterial.friction = 0.6f;
		int kind = slot % 3;
		if (kind == 0)
		{
			b3Sphere s = { { 0.0f, 0.0f, 0.0f }, 0.3f };
			b3CreateSphereShape(body, &sd, &s);
		}
		else if (kind == 1)
		{
			b3Capsule c = { { 0.0f, -0.3f, 0.0f }, { 0.0f, 0.3f, 0.0f }, 0.2f };
			b3CreateCapsuleShape(body, &sd, &c);
		}
		else
		{
			b3BoxHull hull = b3MakeCubeHull(0.3f);
			b3CreateHullShape(body, &sd, &hull.base);
		}
		ring[slot] = body;
	};
	for (int i = 0; i < count; ++i)
	{
		spawn(i);
	}

	double totalMs = 0.0, maxMs = 0.0;
	int cursor = 0;
	const int steps = 240;
	for (int s = 0; s < steps; ++s)
	{
		auto t0 = Clock::now();
		for (int k = 0; k < churnPerStep; ++k)
		{
			b3DestroyBody(ring[cursor]);
			spawn(cursor);
			cursor = (cursor + 1) % count;
		}
		b3World_Step(world, 1.0f / 60.0f, substeps);
		auto t1 = Clock::now();
		double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
		totalMs += ms;
		if (ms > maxMs)
		{
			maxMs = ms;
		}
	}
	printf("churn %d/%d sub=%d       | avg %7.3f ms  max %7.3f ms\n", count, churnPerStep, substeps, totalMs / steps,
		   maxMs);
	b3DestroyWorld(world);
}

int main(int argc, char** argv)
{
	int substeps = argc > 1 ? atoi(argv[1]) : 4;
	printf("box3d benchmark (single thread), substeps=%d\n\n", substeps);
	BenchPyramid(20, substeps);
	BenchStackGrid(8, 6, substeps);
	BenchMixedRain(1000, substeps);
	BenchMiniTowers(20, 3, substeps);
	BenchChurn(1500, 25, substeps);
	return 0;
}
