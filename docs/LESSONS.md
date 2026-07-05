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
