# Field notes: building a stable XPBD engine

Lessons learned the hard way while bringing Melon3D from "towers collapse
after 20 seconds" to "pyramids sleep within millimeters of their spawn
pose". Each of these cost a debugging session; recorded so they cost you
nothing.

## Solver

**Velocity-pass restitution must be unconditional.** Applying the
restitution/normal-velocity fix (paper eq. 29) only when it increases
separation seems safer, but that clamp is exactly what lets the fictitious
separating velocity created by positional projection survive. Symptom:
deep stacks chatter at 0.2–0.7 m/s forever and never sleep.

**Static friction placement decides whether towers live.** Interleaved
with the normal Gauss-Seidel loop, per-point friction makes stacks *walk*
sideways (later normal corrections keep reintroducing tangential drift).
Moved after the normal iterations at full gain, it acts as a base
re-centering controller with loop gain > 1 and *pumps* the tower's sway
mode — towers self-destruct after ~20 simulated seconds. The stable
configuration: under-relaxed corrections (×0.5) toward **persistent
stiction anchors**.

**Stiction anchors must survive manifold regeneration.** If friction
re-anchors every time the manifold is rebuilt, micro-creep ratchets the
target forward and piles never truly stop (they also never sleep, which
dominates benchmark averages). Matching points by feature id and
inheriting the anchors — with a slip release at 2 cm of stretch or a
2× cone overload — arrests creep completely. This single change took the
pyramid scene from "drifts a meter and never sleeps" to "top box within
3 mm of ideal, asleep by step 60".

**Freeze responses, never residuals.** Freezing the effective normal mass
and the angular response Jacobian `J = I⁻¹(r×n)` per substep is safe and
fast. Freezing the *anchors used to measure separation* is not: in a
positional solver the rotational part of a correction only shows up in
the residual through re-rotated anchors, so frozen residual anchors
systematically over-correct until stacks explode.

**First-order quaternion renormalization is enough in the hot path.**
`q · (1.5 − 0.5·|q|²)` after each tiny correction, with one exact
normalize per substep at integration. Error is O(ε²); the sqrt removal is
measurable because corrections run ~6 times per contact point per substep.

**Cap depenetration speed** (3 m/s here). An uncapped positional
projection acts like an impulsive kick when a rocking edge digs in.

## Collision & broadphase

**Per-substep narrowphase wants a pose-delta cache.** Contacts whose
bodies moved less than 1 mm / 1 mrad reuse their manifold. At rest this
makes narrowphase nearly free, and the cached anchors double as the
stiction targets above.

**Moved-only pair updates need a predictive breach test.** With
persistent fat AABBs, re-pairing only bodies that left their fat box is a
big win — but the containment test must use the *one-step-ahead* position
(tight AABB extended by `v·dt`). Testing the current tight box re-pairs
fast bodies one step too late; a 50 m/s projectile then meets a thin wall
with a deep-contact push-out from the wrong side and tunnels.

## Multithreading

**A pure lock-free fork-join pool has a straggler race.** A worker
delayed between reading the job id and joining can observe the *next*
job's reset work cursor combined with the *old* job's bounds and call a
dangling `std::function` out of range (manifested as a rare access
violation under heavy scene churn). Fix: job *discovery* and *joining*
stay lock-free spins, but job *entry* and *retirement* are handshaked
through a mutex-guarded `jobOpen` window.

**Graph coloring pays only above a size threshold.** With ~40 contacts
per color, the per-color dispatch barrier costs more than the work. Keep
the color count low (8 here), send leftovers to a serial overflow pass,
route islands under ~300 contacts to the plain per-island path — and pick
the path by island *size*, never by thread count, so results stay
bit-identical across worker counts.

**Sorted-descending islands must be dispatched at grain 1.** Batching
work items when the list is sorted biggest-first hands all the big
islands to one worker.

**`std::vector<bool>` is not thread-safe even for disjoint indices.**
Parallel writers need `std::vector<uint8_t>`.

## Real-world realism (the v1.2 friction and rotation audit)

**Impulse budgets and velocities are different units.** The dynamic
friction clamp compared `mu * lambdaN * invH` (an impulse) against the
slip *velocity* and then divided by the effective mass again — inflating
kinetic friction by 1/w, two orders of magnitude for anything heavier
than a few kilograms. Every slider was effectively glued: boxes stuck to
slopes far beyond the friction cone and even fell asleep mid-slide.
Symptom-level masking (rolling resistance) had been hiding part of this.
Dimensional analysis on every solver clamp is worth an audit pass.

**A rolling body's stored contact anchor migrates.** The manifold keeps
material-point anchors; on a curved shape the true contact point moves
around the surface as the body rolls, so by freeze time the anchor sits
~w·h off the bottom and `r x n != 0`. Every normal impulse then applies a
phantom torque about the roll axis *proportional to spin*: a frictionless
sphere spinning in place accelerated its own spin at +20%/s, exactly
g·h·m·r/I. Response Jacobians for curved shapes must use the *geometric*
arm (center -> surface along the normal), not the stored material point.

**Positional stiction needs two gates.** It exists to pin resting face
contacts (towers), but it must not touch (a) sliding contacts — gate by
contact tangential speed, or reconstruction turns Coulomb sliding into a
1 cm/s crawl that falls asleep on a 35-degree slope; (b) curved shapes -
their contact point migrates, and the lagging pull-back torque pumps
rollers. Real spheres hold on slopes through rolling resistance, not
tangential stiction.

**Pile noise and genuine sliding are told apart by direction
coherence.** A confined pile recirculates micro-slip through positional
projections faster than the per-substep Coulomb budget can drain it — so
300-body piles never slept. But that noise flips direction chaotically,
while a real slope onset pushes the same way every substep. Storing the
last slip direction per manifold point and applying an asperity-style
full stop only to slow *incoherent* slip freezes piles without touching
Coulomb onset (the effective static-vs-kinetic ratio lands at ~1.1-1.15,
matching real materials where mu_s > mu_k).

**First-order quaternion updates bleed angular momentum.** The
`q + 0.5 h w q`-then-normalize integrator under-rotates by (wh)^2/12 per
substep, and the matching small-angle velocity reconstruction compounds
it: a fast tumble lost 26% of |L| over 10 seconds. An exact
exponential-map update paired with an exact log-map reconstruction cuts
that to ~7% (the remainder is the implicit gyroscopic solve, which
dissipates by design - it never gains energy).

**Cross-build bit-identity is not a valid acceptance gate under FMA.**
With `-O3 -mavx2 -mfma`, adding any code to a hot translation unit
legally changes which expressions the compiler contracts into fused
multiply-adds, which changes rounding, which changes trajectories.
Proof technique: rebuild both versions with `-ffp-contract=off` and
compare those checksums instead; within-build 1T-vs-NT determinism is
unaffected and always valid.
