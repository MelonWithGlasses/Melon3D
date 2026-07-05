// Melon3D benchmarks: pyramid stack, box rain, mixed shapes.
// Prints per-stage timings from the built-in profiler.
#include "melon3d.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <thread>
#include <chrono>
#include <vector>

static void AddGround(m3d_world* world)
{
	m3d_body_def bd = m3d_body_def_default();
	bd.type = M3D_BODY_STATIC;
	bd.position = m3d_v3(0.0f, -1.0f, 0.0f);
	m3d_shape_def sd = m3d_shape_def_default();
	sd.type = M3D_SHAPE_BOX;
	sd.halfExtents = m3d_v3(200.0f, 1.0f, 200.0f);
	m3d_body_create(world, &bd, &sd);
}

static void RunScene(const char* name, m3d_world* world, int steps, int substeps)
{
	float totalMs = 0.0f;
	float maxMs = 0.0f;
	float broadMs = 0.0f, narrowMs = 0.0f, solveMs = 0.0f;
	for (int i = 0; i < steps; ++i)
	{
		m3d_world_step(world, 1.0f / 60.0f, substeps);
		m3d_profile p = m3d_world_profile(world);
		totalMs += p.stepMs;
		broadMs += p.broadphaseMs;
		narrowMs += p.narrowphaseMs;
		solveMs += p.solverMs;
		if (p.stepMs > maxMs)
		{
			maxMs = p.stepMs;
		}
	}
	m3d_profile p = m3d_world_profile(world);
	printf("%-24s bodies=%5d contacts=%5d islands=%4d | avg %7.3f ms (bp %6.3f, np %6.3f, solve %6.3f) max %7.3f ms\n",
		   name, p.bodyCount, p.touchingContactCount, p.islandCount, totalMs / (float)steps, broadMs / (float)steps,
		   narrowMs / (float)steps, solveMs / (float)steps, maxMs);
}

// Classic box pyramid: base*base/2 boxes leaning on each other in one island.
static void BenchPyramid(int workerCount, int base, int substeps)
{
	m3d_world_def wd = m3d_world_def_default();
	wd.workerCount = workerCount;
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	const float he = 0.5f;
	for (int row = 0; row < base; ++row)
	{
		int countInRow = base - row;
		float y = he + 2.0f * he * (float)row;
		float x0 = -he * (float)(countInRow - 1);
		for (int i = 0; i < countInRow; ++i)
		{
			m3d_body_def bd = m3d_body_def_default();
			bd.type = M3D_BODY_DYNAMIC;
			bd.position = m3d_v3(x0 + 2.0f * he * (float)i, y, 0.0f);
			m3d_shape_def sd = m3d_shape_def_default();
			sd.type = M3D_SHAPE_BOX;
			sd.halfExtents = m3d_v3(he, he, he);
			sd.friction = 0.6f;
			m3d_body_create(world, &bd, &sd);
		}
	}

	char name[64];
	snprintf(name, sizeof(name), "pyramid(%d) sub=%d %dT", base, substeps, workerCount);
	RunScene(name, world, 240, substeps);
	m3d_world_destroy(world);
}

// Grid of separated small stacks: many independent islands (parallel-friendly).
static void BenchStackGrid(int workerCount, int grid, int height, int substeps)
{
	m3d_world_def wd = m3d_world_def_default();
	wd.workerCount = workerCount;
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	for (int gx = 0; gx < grid; ++gx)
	{
		for (int gz = 0; gz < grid; ++gz)
		{
			float x = 4.0f * (float)gx - 2.0f * (float)grid;
			float z = 4.0f * (float)gz - 2.0f * (float)grid;
			for (int j = 0; j < height; ++j)
			{
				m3d_body_def bd = m3d_body_def_default();
				bd.type = M3D_BODY_DYNAMIC;
				bd.position = m3d_v3(x, 0.5f + 1.0f * (float)j, z);
				m3d_shape_def sd = m3d_shape_def_default();
				sd.type = M3D_SHAPE_BOX;
				sd.halfExtents = m3d_v3(0.5f, 0.5f, 0.5f);
				m3d_body_create(world, &bd, &sd);
			}
		}
	}

	char name[64];
	snprintf(name, sizeof(name), "stacks %dx%dx%d sub=%d %dT", grid, grid, height, substeps, workerCount);
	RunScene(name, world, 240, substeps);
	m3d_world_destroy(world);
}

// Mixed shape rain onto the ground.
static void BenchMixedRain(int workerCount, int count, int substeps)
{
	m3d_world_def wd = m3d_world_def_default();
	wd.workerCount = workerCount;
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	unsigned seed = 12345;
	auto frand = [&seed]() {
		seed = seed * 1664525u + 1013904223u;
		return (float)(seed >> 8) / (float)0x00FFFFFF;
	};

	for (int i = 0; i < count; ++i)
	{
		m3d_body_def bd = m3d_body_def_default();
		bd.type = M3D_BODY_DYNAMIC;
		bd.position = m3d_v3(20.0f * frand() - 10.0f, 1.0f + 0.1f * (float)i, 20.0f * frand() - 10.0f);
		m3d_shape_def sd = m3d_shape_def_default();
		int kind = i % 3;
		if (kind == 0)
		{
			sd.type = M3D_SHAPE_SPHERE;
			sd.radius = 0.3f;
		}
		else if (kind == 1)
		{
			sd.type = M3D_SHAPE_CAPSULE;
			sd.radius = 0.2f;
			sd.halfHeight = 0.3f;
		}
		else
		{
			sd.type = M3D_SHAPE_BOX;
			sd.halfExtents = m3d_v3(0.3f, 0.3f, 0.3f);
		}
		m3d_body_create(world, &bd, &sd);
	}

	char name[64];
	snprintf(name, sizeof(name), "rain %d sub=%d %dT", count, substeps, workerCount);
	RunScene(name, world, 300, substeps);
	m3d_world_destroy(world);
}

// Field of small independent towers: many islands, the architecture's
// best case for island-parallel solving.
static void BenchMiniTowers(int workerCount, int grid, int height, int substeps)
{
	m3d_world_def wd = m3d_world_def_default();
	wd.workerCount = workerCount;
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	for (int gx = 0; gx < grid; ++gx)
	{
		for (int gz = 0; gz < grid; ++gz)
		{
			float x = 3.0f * (float)gx - 1.5f * (float)grid;
			float z = 3.0f * (float)gz - 1.5f * (float)grid;
			for (int j = 0; j < height; ++j)
			{
				m3d_body_def bd = m3d_body_def_default();
				bd.type = M3D_BODY_DYNAMIC;
				bd.position = m3d_v3(x, 0.5f + 1.0f * (float)j, z);
				m3d_shape_def sd = m3d_shape_def_default();
				sd.type = M3D_SHAPE_BOX;
				sd.halfExtents = m3d_v3(0.5f, 0.5f, 0.5f);
				m3d_body_create(world, &bd, &sd);
			}
		}
	}

	char name[64];
	snprintf(name, sizeof(name), "towers %dx%dx%d sub=%d %dT", grid, grid, height, substeps, workerCount);
	RunScene(name, world, 240, substeps);
	m3d_world_destroy(world);
}

// Streaming churn: bodies are continuously destroyed and respawned while
// raining down - debris/projectile workload. Create/destroy time is part
// of the measurement (timed externally around the whole frame).
static void BenchChurn(int workerCount, int count, int churnPerStep, int substeps)
{
	m3d_world_def wd = m3d_world_def_default();
	wd.workerCount = workerCount;
	m3d_world* world = m3d_world_create(&wd);
	AddGround(world);

	unsigned seed = 777;
	auto frand = [&seed]() {
		seed = seed * 1664525u + 1013904223u;
		return (float)(seed >> 8) / (float)0x00FFFFFF;
	};

	std::vector<m3d_body_id> ring((size_t)count);
	auto spawn = [&](int slot) {
		m3d_body_def bd = m3d_body_def_default();
		bd.type = M3D_BODY_DYNAMIC;
		bd.position = m3d_v3(24.0f * frand() - 12.0f, 6.0f + 8.0f * frand(), 24.0f * frand() - 12.0f);
		m3d_shape_def sd = m3d_shape_def_default();
		int kind = slot % 3;
		if (kind == 0)
		{
			sd.type = M3D_SHAPE_SPHERE;
			sd.radius = 0.3f;
		}
		else if (kind == 1)
		{
			sd.type = M3D_SHAPE_CAPSULE;
			sd.radius = 0.2f;
			sd.halfHeight = 0.3f;
		}
		else
		{
			sd.type = M3D_SHAPE_BOX;
			sd.halfExtents = m3d_v3(0.3f, 0.3f, 0.3f);
		}
		ring[slot] = m3d_body_create(world, &bd, &sd);
	};
	for (int i = 0; i < count; ++i)
	{
		spawn(i);
	}

	using Clock = std::chrono::steady_clock;
	double totalMs = 0.0, maxMs = 0.0;
	int cursor = 0;
	const int steps = 240;
	for (int s = 0; s < steps; ++s)
	{
		auto t0 = Clock::now();
		for (int k = 0; k < churnPerStep; ++k)
		{
			m3d_body_destroy(world, ring[cursor]);
			spawn(cursor);
			cursor = (cursor + 1) % count;
		}
		m3d_world_step(world, 1.0f / 60.0f, substeps);
		auto t1 = Clock::now();
		double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
		totalMs += ms;
		if (ms > maxMs)
		{
			maxMs = ms;
		}
	}
	printf("churn %d/%d sub=%d %dT%*s | avg %7.3f ms  max %7.3f ms\n", count, churnPerStep, substeps, workerCount,
		   workerCount > 9 ? 0 : 1, "", totalMs / steps, maxMs);
	m3d_world_destroy(world);
}

int main(int argc, char** argv)
{
	int hw = (int)std::thread::hardware_concurrency();
	if (hw < 1)
	{
		hw = 4;
	}
	int threads = argc > 1 ? atoi(argv[1]) : hw;
	printf("Melon3D benchmark, hardware threads: %d, using up to %d\n\n", hw, threads);

	for (int sub = 4; sub <= 8; sub += 4)
	{
		BenchPyramid(1, 20, sub);
		if (threads > 1)
		{
			BenchPyramid(threads, 20, sub);
		}
		BenchStackGrid(1, 8, 6, sub);
		if (threads > 1)
		{
			BenchStackGrid(threads, 8, 6, sub);
		}
		BenchMixedRain(1, 1000, sub);
		if (threads > 1)
		{
			BenchMixedRain(threads, 1000, sub);
		}
		if (sub == 4)
		{
			BenchMiniTowers(1, 20, 3, sub);
			if (threads > 1)
			{
				BenchMiniTowers(threads, 20, 3, sub);
			}
			BenchChurn(1, 1500, 25, sub);
			if (threads > 1)
			{
				BenchChurn(threads, 1500, 25, sub);
			}
		}
	}

	return 0;
}
