# Melon3D developer guide

Practical recipes for everything the engine does. All snippets are plain C
(the API is C89-compatible; the implementation is C++17).

- [Integrating the library](#integrating-the-library)
- [World and stepping](#world-and-stepping)
- [Bodies and shapes](#bodies-and-shapes)
- [Materials](#materials)
- [Collision filtering](#collision-filtering)
- [Sensors (trigger volumes)](#sensors-trigger-volumes)
- [Forces, impulses, velocities](#forces-impulses-velocities)
- [Joints](#joints)
- [Queries: rays, sphere casts, overlaps](#queries-rays-sphere-casts-overlaps)
- [Contact events](#contact-events)
- [A minimal character controller](#a-minimal-character-controller)
- [Threading and determinism](#threading-and-determinism)
- [Sleeping](#sleeping)
- [Profiling and tuning](#profiling-and-tuning)

## Integrating the library

Three ways, pick one:

**CMake subdirectory / FetchContent** (recommended):

```cmake
include(FetchContent)
FetchContent_Declare(melon3d
    GIT_REPOSITORY https://github.com/MelonWithGlasses/Melon3D.git
    GIT_TAG main)
FetchContent_MakeAvailable(melon3d)
target_link_libraries(my_game PRIVATE melon3d::melon3d)
```

Tests and benchmarks are skipped automatically when Melon3D is not the
top-level project.

**Copy the sources.** The whole engine is `include/melon3d.h` plus the four
files in `src/`. Add them to any build system, compile as C++17, done. No
dependencies beyond the C++ standard library and threads.

**Single file.** `python tools/amalgamate.py` writes
`dist/melon3d_single.cpp`; drop it and `melon3d.h` into your project and
compile the one `.cpp`. Good for jam games and quick experiments.

On x86, build with AVX2 for the SIMD contact packets (`-mavx2 -mfma` on
GCC/Clang, `/arch:AVX2` on MSVC). Everything works without it, just slower
on big piles. ARM builds autovectorize with NEON out of the box.

## World and stepping

```c
m3d_world_def wd = m3d_world_def_default();
wd.gravity = m3d_v3(0.0f, -9.81f, 0.0f);
wd.workerCount = 4;          // worker threads incl. the caller; results are
                             // bit-identical for ANY value here
m3d_world* world = m3d_world_create(&wd);

// fixed timestep, 4 substeps is the sweet spot for stacks
m3d_world_step(world, 1.0f / 60.0f, 4);

m3d_world_destroy(world);
```

Use a fixed `dt`. XPBD substeps replace solver iterations: more substeps =
stiffer stacks and joints (8 for tall towers, 2 for loose debris is fine).
`m3d_world_set_gravity` changes gravity at runtime (sleeping bodies stay
asleep until something wakes them).

## Bodies and shapes

One shape per body: sphere, capsule (Y axis), or box.

```c
m3d_body_def bd = m3d_body_def_default();    // static by default
bd.type = M3D_BODY_DYNAMIC;                  // or M3D_BODY_KINEMATIC
bd.position = m3d_v3(0, 5, 0);
bd.orientation = m3d_quat_axis_angle(m3d_v3(0, 0, 1), 0.3f);
bd.userData = (uint64_t)(uintptr_t)my_entity; // round-trips untouched

m3d_shape_def sd = m3d_shape_def_default();
sd.type = M3D_SHAPE_CAPSULE;
sd.radius = 0.4f;
sd.halfHeight = 0.5f;      // inner segment half-length; total height = 1.8
sd.density = 1000.0f;      // kg/m^3 - mass comes from the shape volume

m3d_body_id id = m3d_body_create(world, &bd, &sd);
```

Ids are generational: after `m3d_body_destroy` the id goes stale and every
accessor returns a safe default (`m3d_body_is_valid` tells you explicitly).
Kinematic bodies follow their velocities exactly, push dynamic bodies, and
are never pushed back — set velocities, don't set transforms every frame.
`m3d_body_set_transform` is a teleport: fine for spawning and cutscenes.

## Materials

```c
sd.friction = 0.6f;           // Coulomb, combined as sqrt(a*b)
sd.restitution = 0.3f;        // bounciness, combined as max(a,b)
sd.rollingResistance = 0.02f; // spheres/capsules stop rolling; 0 = billiards
```

Runtime changes (`m3d_body_set_friction`, `m3d_body_set_restitution`) also
update already-existing contacts, so a "now the floor is ice" switch works
immediately.

## Collision filtering

Box2D-compatible semantics:

```c
enum { LAYER_WORLD = 1 << 0, LAYER_PLAYER = 1 << 1, LAYER_ENEMY = 1 << 2,
       LAYER_DEBRIS = 1 << 3 };

// debris collides with the world but not with players or other debris
sd.filter.categoryBits = LAYER_DEBRIS;
sd.filter.maskBits = LAYER_WORLD;

// ragdoll parts of one character: same negative group = never collide
sd.filter.groupIndex = -playerId;
```

Two bodies collide when they share a positive `groupIndex` (always), do NOT
share a negative one (never), or otherwise when
`(a.category & b.mask) && (b.category & a.mask)`. Filtering happens at pair
creation, so filtered pairs cost nothing — no manifold, no memory.

Bodies connected by a joint never collide with each other.

## Sensors (trigger volumes)

```c
m3d_shape_def zone = m3d_shape_def_default();
zone.type = M3D_SHAPE_BOX;
zone.halfExtents = m3d_v3(2, 1, 2);
zone.isSensor = true;
m3d_body_id trigger = m3d_body_create(world, &zd, &zone);  // static body
```

Sensors detect overlap and emit contact begin/end events (flagged
`isSensor`), but generate no forces, never wake sleeping bodies, and are
excluded from the solver entirely. Rays and sphere casts *do* hit sensors —
give them their own category and mask them out of the cast when needed.

## Forces, impulses, velocities

```c
// forces accumulate and act for exactly the next step - call every frame
// for a sustained effect (hover, wind, thrust)
m3d_body_apply_force(world, id, m3d_v3(0, 900.0f, 0));
m3d_body_apply_force_at_point(world, id, thrust, nozzleWorldPos);
m3d_body_apply_torque(world, id, m3d_v3(0, 50.0f, 0));

// impulses change velocity immediately (explosions, jumps, hits)
m3d_body_apply_impulse(world, id, m3d_v3(0, 5.0f * mass, 0), hitPoint);
m3d_body_apply_angular_impulse(world, id, m3d_v3(0, 2.0f, 0));

// or set velocities outright (character controllers do this)
m3d_body_set_linear_velocity(world, id, m3d_v3(move.x, vy, move.z));
```

All of these wake the body. `bd.gravityScale = 0` makes a body ignore
gravity (useful for projectiles with custom ballistics).

## Joints

All joints are solved positionally (XPBD), so they do not drift apart under
load; more substeps = stiffer.

```c
// distance: rigid rod, or a spring via hertz/dampingRatio
m3d_distance_joint_def dd = { bodyA, bodyB, anchorInA, anchorInB,
                              /*length*/ 0.0f,   // 0 = current distance
                              /*hertz*/ 4.0f, /*dampingRatio*/ 0.5f };
m3d_joint_create_distance(world, &dd);

// ball: anchors coincide, rotation free (ragdoll shoulders)
m3d_ball_joint_def bj = { bodyA, bodyB, anchorInA, anchorInB };
m3d_joint_create_ball(world, &bj);

// hinge: one rotational DOF, optional limit and motor (doors, wheels)
m3d_hinge_joint_def hd = m3d_hinge_joint_def_default();
hd.bodyA = doorFrame;
hd.bodyB = door;
hd.localAnchorA = m3d_v3(0.5f, 0, 0);   // pivot on the frame
hd.localAnchorB = m3d_v3(-0.5f, 0, 0);  // hinge edge of the door
hd.localAxisA = m3d_v3(0, 1, 0);        // swing around vertical
// localAxisB stays zero -> derived from the current poses
hd.enableLimit = true;
hd.lowerAngle = -0.5f * M3D_PI;
hd.upperAngle = 0.0f;
m3d_joint_id hinge = m3d_joint_create_hinge(world, &hd);

// wheels: no limit, motor on
m3d_joint_hinge_set_motor(world, hinge, true, /*rad/s*/ 20.0f, /*N*m*/ 150.0f);
float angle = m3d_joint_hinge_angle(world, hinge);  // radians, 0 at creation

// weld: full lock (breakable structures - destroy the joint to "break" it)
m3d_weld_joint_def wj = m3d_weld_joint_def_default();
wj.bodyA = a; wj.bodyB = b;
wj.localAnchorA = m3d_v3(0.5f, 0, 0);
wj.localAnchorB = m3d_v3(-0.5f, 0, 0);
m3d_joint_create_weld(world, &wj);
```

The hinge motor is torque-limited: with `motorSpeed = 0` it acts as a
friction brake. Angles are measured around the hinge axis, zero at the pose
the joint was created in.

## Queries: rays, sphere casts, overlaps

```c
// closest-hit ray
m3d_ray_result r = m3d_world_ray_cast(world, muzzle, m3d_scale(100.0f, dir));
if (r.hit) { /* r.body, r.point, r.normal, r.fraction */ }

// same but only against chosen categories
r = m3d_world_ray_cast_filtered(world, muzzle, dir100, LAYER_WORLD | LAYER_ENEMY);

// swept sphere - the workhorse of character movement
r = m3d_world_sphere_cast(world, feet, 0.4f, m3d_v3(0, -0.6f, 0), LAYER_WORLD);

// all bodies in a region (explosion radius, selection box)
bool visit(m3d_body_id id, void* ctx) { /* ... */ return true; }
m3d_aabb box = { m3d_v3(-5, 0, -5), m3d_v3(5, 10, 5) };
m3d_world_overlap_aabb(world, box, visit, &myCtx);
```

Sphere casts are exact against spheres and capsules and use conservative
advancement against boxes. Returning `false` from the overlap callback
stops the query.

## Contact events

```c
m3d_contact_events ev = m3d_world_contact_events(world);
for (int i = 0; i < ev.beginCount; ++i) {
    const m3d_contact_event* e = &ev.beginEvents[i];
    if (e->isSensor) { on_trigger_enter(e->bodyA, e->bodyB); }
    else            { on_touch(e->bodyA, e->bodyB); }
}
```

Events are valid until the next `m3d_world_step`. Use `userData` to map
bodies back to your entities.

## A minimal character controller

A kinematic capsule moved by sphere casts — robust and completely under
your control:

```c
// probe the ground
m3d_ray_result ground = m3d_world_sphere_cast(world, pos, capRadius,
                                              m3d_v3(0, -stepDown, 0), LAYER_WORLD);
bool grounded = ground.hit && ground.normal.y > 0.7f;

// horizontal move, sliding on hit
m3d_vec3 delta = m3d_scale(dt, wishVelocity);
m3d_ray_result wall = m3d_world_sphere_cast(world, pos, capRadius, delta, LAYER_WORLD);
if (wall.hit) {
    // remove the blocked component, keep the tangential one
    m3d_vec3 n = wall.normal;
    delta = m3d_sub(delta, m3d_scale(m3d_dot(delta, n), n));
}
m3d_body_set_linear_velocity(world, character, m3d_scale(1.0f / dt, delta));
```

Give the character its own filter category so casts can exclude it, and
keep it kinematic so piles of debris cannot shove it around (or dynamic if
you want that).

## Threading and determinism

`wd.workerCount = N` uses the built-in fork-join pool. **Simulation results
are bit-identical for any worker count** — parallel paths are chosen by
workload shape (island size), never by thread count. The same binary on the
same machine with the same inputs replays exactly; lockstep networking over
identical builds works.

The API is not thread-safe: call all `m3d_*` functions for one world from
one thread (the pool parallelizes inside the step on its own).

## Sleeping

Islands that stay quiet for 0.5 s fall asleep and cost almost nothing —
sleeping bodies migrate to a *stable* broadphase tier that is not even
scanned. Tens of thousands of resting bodies are fine. Anything that
touches a sleeping body (contact from an awake body, joint, force, setter)
wakes it, and wake propagates through the island. `bd.enableSleep = false`
opts a body out.

## Profiling and tuning

```c
m3d_profile p = m3d_world_profile(world);
// p.stepMs, broadphaseMs, narrowphaseMs, solverMs,
// p.bodyCount, awakeBodyCount, contactCount, islandCount
```

Rules of thumb:

- `wd.gridCellSize` ≈ the diameter of your median object (default 1.5 m).
- Substeps: 4 default, 8 for tall towers/strong motors, 2 for debris.
- Prefer several medium boxes over one extremely long thin one.
- Give bullets a small `radius` and high velocity — speculative margins
  handle continuous collision automatically (see the no-tunneling test).
