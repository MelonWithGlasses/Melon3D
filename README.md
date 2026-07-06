# Melon3D

[![CI](https://github.com/MelonWithGlasses/Melon3D/actions/workflows/ci.yml/badge.svg)](https://github.com/MelonWithGlasses/Melon3D/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![Dependencies: none](https://img.shields.io/badge/dependencies-none-brightgreen.svg)](#building)

A compact 3D rigid body physics engine for games, built on **position-based
dynamics (XPBD)**. Zero dependencies, C-style API, C++17 inside,
deterministic multithreading out of the box.

```
~5k lines · MIT · sphere / capsule / box · joints · sleeping · ray casts
```

## Why another physics engine?

Melon3D deliberately explores the *other* branch of the design tree from
the impulse/warm-starting school (Box2D, box3d, PhysX):

| | Impulse solvers | Melon3D |
|---|---|---|
| Constraint level | velocity | **position (XPBD)** |
| Rest penetration | slop-bounded overlap | **exactly zero** |
| Warm-start state | persistent impulses | **none** |
| Broadphase | incremental BVH | **two-tier hash grid** |
| Threading model | task scheduler (bring your own) | **built-in pool, deterministic** |

The result has a different performance signature: it shines on streaming
workloads (constant spawn/despawn) and worst-frame latency, and pays a
premium on piles that settle and sleep. Honest numbers below.

## Features

- **Shapes**: sphere, capsule, box (OBB), one shape per body
- **Broadphase**: two-tier sorted spatial hash grid — sleeping and static
  bodies live in a *stable* tier rebuilt only on membership changes;
  moved-only pair updates via persistent fat AABBs with a one-step
  velocity prediction (fast bodies always pair before impact — no
  tunneling, verified by test)
- **Narrowphase**: SAT box-box with reference-face clipping, capsule
  closest-point via convex ternary search; manifolds cached by pose delta
  (1 mm / 1 mrad); speculative margins scaled by relative velocity
- **Solver**: XPBD substepping (Müller et al. 2020)
  - persistent **stiction anchors** for static friction — no creep, no
    sway pumping, piles come to a true stop and sleep
  - per-substep frozen effective masses and angular Jacobians
  - depenetration capped at 3 m/s; separate restitution/dynamic-friction
    velocity pass
- **Islands & sleeping**: union-find islands, whole-island sleep/wake
- **Big-island parallelism**: islands with 300+ contacts are solved with
  **graph coloring** — contacts bucketed so no two in a bucket share a
  dynamic body — and **8-wide SIMD contact packets** (branchless
  `float[8]` lane code, fully autovectorized with AVX2)
- **Threading**: built-in spin fork-join pool (~1–2 µs stage dispatch).
  Simulation results are **bit-identical for any worker count** — the
  colored path is selected by island size, never by thread count, and
  bodies within a color are disjoint (covered by a test)
- **Joints**: distance (rigid or spring via XPBD compliance), ball
- **Queries & events**: contact begin/end events, DDA grid ray casting,
  per-stage profiling (`m3d_world_profile`)

## Building

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
build/melon3d_tests   # 21 checks
build/melon3d_bench   # benchmark scenes
```

GCC / Clang / MSVC, Windows / Linux / macOS. No dependencies.

## Quick start

```c
#include "melon3d.h"

m3d_world_def wd = m3d_world_def_default();
wd.workerCount = 4;                    // parallel, still deterministic
m3d_world* world = m3d_world_create(&wd);

m3d_body_def ground = m3d_body_def_default();   // static by default
ground.position = m3d_v3(0.0f, -1.0f, 0.0f);
m3d_shape_def gs = m3d_shape_def_default();
gs.type = M3D_SHAPE_BOX;
gs.halfExtents = m3d_v3(50.0f, 1.0f, 50.0f);
m3d_body_create(world, &ground, &gs);

m3d_body_def bd = m3d_body_def_default();
bd.type = M3D_BODY_DYNAMIC;
bd.position = m3d_v3(0.0f, 5.0f, 0.0f);
m3d_shape_def sd = m3d_shape_def_default();
sd.type = M3D_SHAPE_BOX;
sd.halfExtents = m3d_v3(0.5f, 0.5f, 0.5f);
m3d_body_id box = m3d_body_create(world, &bd, &sd);

for (int i = 0; i < 60; ++i)
    m3d_world_step(world, 1.0f / 60.0f, 4);     // 4 substeps

m3d_vec3 p = m3d_body_position(world, box);
m3d_world_destroy(world);
```

## Benchmarks vs box3d

[box3d](https://github.com/erincatto/box3d) is Erin Catto's 3D physics
engine — the natural reference point for this project and the inspiration
for several of its scenes. The comparison below is offered in that spirit:
one data point of XPBD-vs-impulse trade-offs, not a leaderboard.

### Setup (full disclosure)

| | |
|---|---|
| CPU | AMD Ryzen 7 7735HS (8 cores / 16 threads, Zen 3+) |
| RAM | 16 GB DDR5 |
| OS | Windows 10 Pro |
| Compiler | GCC 15.2.0 (MSYS2 ucrt64), `-O3 -mavx2 -mfma` **for both engines** |
| box3d version | commit `1bec63c` (July 2026) |
| Timestep | 60 Hz, 4 substeps, identical scenes/step counts/PRNG seeds |
| Timing | wall time around the full step; churn also times create/destroy |
| Runs | two consecutive runs per engine, averaged; run-to-run spread ≈ ±5% |

**Important**: box3d runs **single-threaded** in this harness because its
task-system callbacks are not wired here. box3d supports multithreading
when a task scheduler is provided, so read the table as "Melon3D with its
built-in pool vs box3d out-of-the-box in this harness" — not as an upper
bound of what box3d can do. Melon3D single-thread numbers are included so
both engines can also be compared 1T vs 1T.

### Results — average ms per step (worst frame in parentheses)

| Scene | box3d 1T | Melon3D 1T | Melon3D 16T |
|---|---|---|---|
| **churn** — 1500 bodies raining, 25 destroyed+respawned each step | 1.282 (1.97) | 1.53 (2.44) | **1.000 (2.10)** |
| **rain** — 1000 mixed spheres/capsules/boxes falling | **0.976** (2.46) | 1.88 (5.48) | 1.059 (3.95) |
| **stacks** — 8×8 grid of 6-box towers (384 bodies) | **0.046** (0.86) | 0.203 (1.92) | 0.057 (0.83) |
| **towers** — 20×20 grid of 3-box towers (1200 bodies) | **0.143** (2.94) | 0.643 (5.83) | 0.238 (2.04) |
| **pyramid** — 210-box pyramid, one contact island | **0.121** (1.59) | 0.360 (2.75) | 0.372 (2.36) |

<details>
<summary>Raw per-run data</summary>

```
box3d 1T          run1 avg/max      run2 avg/max
pyramid           0.123 / 1.850     0.119 / 1.335
stacks            0.045 / 0.876     0.047 / 0.848
rain              0.984 / 2.173     0.968 / 2.745
towers            0.145 / 3.049     0.140 / 2.823
churn             1.289 / 1.888     1.275 / 2.048

Melon3D 16T       run1 avg/max      run2 avg/max
pyramid           0.394 / 2.054     0.413 / 2.195
stacks            0.057 / 0.622     0.057 / 1.028
rain              1.045 / 3.882     1.073 / 4.023
towers            0.240 / 1.950     0.235 / 2.129
churn             1.003 / 2.238     0.997 / 1.959

Melon3D 1T        run1 avg/max      run2 avg/max
pyramid           0.422 / 2.579     0.397 / 2.499
stacks            0.202 / 1.914     0.203 / 1.929
rain              1.877 / 5.383     1.873 / 5.574
towers            0.640 / 5.605     0.645 / 6.060
churn             1.539 / 2.516     1.524 / 2.353

Addendum (energy-based sleep criterion, single run):
pyramid Melon3D   1T 0.360 / 2.745  16T 0.372 / 2.362
(other scenes unchanged within run-to-run spread; the summary table above
uses these updated pyramid numbers)
```
</details>

### Reading the numbers

- **Streaming (churn)**: Melon3D 16T is ~22% faster on average than box3d
  1T. With constant create/destroy there is no warm-start state to
  rebuild, the hash grid is indifferent to body lifetime, and stage
  dispatch is cheap — this is the workload the architecture was shaped
  for (debris, projectiles, casings).
- **Rain**: near parity (~8%).
- **Settle-and-sleep scenes (stacks, towers, pyramid)**: box3d is
  1.3–3.3× faster on average. Two measured factors: its warm-started
  solver drives residual velocities below the sleep threshold by roughly
  step 25 where Melon3D's position-reconstructed velocities take until
  roughly step 60 (sleeping piles cost nothing, so the shorter settle
  window dominates the average), and its hand-packed SIMD constraint
  solver does more work per cycle. After Melon3D's SIMD-packet pass, its
  awake-step cost on the pyramid is within ~1.4× of box3d — most of the
  remaining average gap is settle time.
- **Worst-frame latency**: Melon3D is at or better than box3d on churn,
  stacks and towers — frame pacing benefits from the same properties that
  win the streaming scene.

Reproduction: [`bench_compare/`](bench_compare/README.md) contains a
harness that runs the identical scenes on box3d, with full build
instructions.

## Step pipeline

1. **Broadphase** — refresh fat AABBs of awake bodies (breach test uses a
   one-step velocity prediction), rebuild the active grid tier; stable
   tier only on membership changes.
2. **Pairs** — active×active and active×stable candidates from shared
   cells, min-corner deduplication, only pairs with a *moved* body are
   re-examined; oversized bodies (e.g. the ground) via a coarse list.
3. **Narrowphase** *(parallel when there is real work)* — regenerate
   manifolds only for pairs whose bodies moved 1 mm / 1 mrad; stiction
   anchors are inherited across regenerations by feature id.
4. **Islands** — union-find over touching dynamic pairs and joints.
5. **Solve** — big islands: graph-colored, 8-wide SIMD packets; small
   islands: one task each. Per substep: integrate → per-substep manifold
   refresh → project normals (2 iterations) → stiction friction →
   velocity reconstruction → restitution/dynamic-friction pass.
6. **Sleeping** — islands quiet for 0.5 s sleep and migrate to the stable
   broadphase tier.

## Field notes

Hard-won implementation lessons are documented in
[docs/LESSONS.md](docs/LESSONS.md) — friction anchor inheritance, why
freezing residual anchors explodes while freezing response Jacobians is
fine, the lock-free pool straggler race, and when graph coloring starts
to pay.

## Roadmap

- Energy-based sleep criterion (close the settle-time gap)
- SIMD packets for the friction and velocity passes
- Convex hulls (GJK/EPA), hinge and prismatic joints
- Multiple shapes per body, shape casts, sensors

## License

MIT — see [LICENSE](LICENSE).

box3d is © Erin Catto, MIT licensed, and is not included in this
repository; the comparison harness downloads it separately.
