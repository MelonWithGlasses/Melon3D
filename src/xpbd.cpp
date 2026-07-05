// XPBD island solver (Mueller et al., "Detailed Rigid Body Simulation with
// XPBD", 2020). Per substep:
//   1. integrate positions from velocities
//   2. project position constraints (contacts with static friction, joints)
//   3. reconstruct velocities from the position delta
//   4. velocity pass: dynamic friction and restitution
// There are no persistent impulses and no warm starting: penetration is
// removed by positional projection, so stacks rest exactly on the surface.
#include "m3d_internal.h"

namespace m3d
{

// First-order (Newton) renormalization: corrections are tiny, the error is
// O(eps^2), and a full normalize still runs once per substep at integration.
// Avoids a sqrt+div per correction in the hottest loop of the engine.
static inline m3d_quat QuatRenormFast(m3d_quat q)
{
	float s = 1.5f - 0.5f * (q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
	return { q.x * s, q.y * s, q.z * s, q.w * s };
}

// positional correction: impulse-like vector P applied at offset r
static inline void ApplyPositionCorrection(Body& body, m3d_vec3 p, m3d_vec3 r)
{
	if (body.invMass == 0.0f)
	{
		return;
	}
	body.position = m3d_mul_add(body.position, body.invMass, p);
	m3d_vec3 dw = m3d_mat3_mulv(&body.invInertiaWorld, m3d_cross(r, p));
	m3d_quat wq = { dw.x, dw.y, dw.z, 0.0f };
	m3d_quat dq = m3d_quat_mul(wq, body.rotation);
	m3d_quat q = { body.rotation.x + 0.5f * dq.x, body.rotation.y + 0.5f * dq.y, body.rotation.z + 0.5f * dq.z,
				   body.rotation.w + 0.5f * dq.w };
	body.rotation = QuatRenormFast(q);
}

// positional correction of signed magnitude along the contact normal, using
// the per-substep frozen angular Jacobian (no inertia matmul, no cross)
static inline void ApplyAxisCorrection(Body& body, float dLambdaSigned, m3d_vec3 n, m3d_vec3 angJ)
{
	if (body.invMass == 0.0f)
	{
		return;
	}
	body.position = m3d_mul_add(body.position, dLambdaSigned * body.invMass, n);
	m3d_vec3 dw = m3d_scale(dLambdaSigned, angJ);
	m3d_quat wq = { dw.x, dw.y, dw.z, 0.0f };
	m3d_quat dq = m3d_quat_mul(wq, body.rotation);
	m3d_quat q = { body.rotation.x + 0.5f * dq.x, body.rotation.y + 0.5f * dq.y, body.rotation.z + 0.5f * dq.z,
				   body.rotation.w + 0.5f * dq.w };
	body.rotation = QuatRenormFast(q);
}

// velocity impulse of signed magnitude along the normal via the frozen Jacobian
static inline void ApplyAxisVelocity(Body& body, float jSigned, m3d_vec3 n, m3d_vec3 angJ)
{
	if (body.invMass == 0.0f)
	{
		return;
	}
	body.linearVelocity = m3d_mul_add(body.linearVelocity, jSigned * body.invMass, n);
	body.angularVelocity = m3d_mul_add(body.angularVelocity, jSigned, angJ);
}

static inline void ApplyVelocityImpulse(Body& body, m3d_vec3 p, m3d_vec3 r)
{
	// static/kinematic bodies are shared read-only across color buckets:
	// they must not be written, even with a zero delta
	if (body.invMass == 0.0f)
	{
		return;
	}
	body.linearVelocity = m3d_mul_add(body.linearVelocity, body.invMass, p);
	m3d_vec3 dw = m3d_mat3_mulv(&body.invInertiaWorld, m3d_cross(r, p));
	body.angularVelocity = m3d_add(body.angularVelocity, dw);
}

// generalized inverse mass of a correction along `axis` applied at offset r
static inline float GeneralizedInvMass(const Body& b, m3d_vec3 r, m3d_vec3 axis)
{
	m3d_vec3 rn = m3d_cross(r, axis);
	return b.invMass + m3d_dot(rn, m3d_mat3_mulv(&b.invInertiaWorld, rn));
}

// ---------------------------------------------------------------------------
// Contact position projection
// ---------------------------------------------------------------------------

static void SolveContactNormals(m3d_world* world, Contact& c, float maxCorrection)
{
	Body& a = world->bodies[c.bodyA];
	Body& b = world->bodies[c.bodyB];
	Manifold& man = c.manifold;
	m3d_vec3 n = man.normal;

	for (int32_t i = 0; i < man.pointCount; ++i)
	{
		ManifoldPoint& mp = man.points[i];
		if (mp.wN <= 0.0f)
		{
			continue;
		}

		// anchors must be re-rotated every visit: the rotational part of a
		// projection only shows up in the residual through fresh anchors
		// (frozen anchors systematically over-correct and explode stacks);
		// the effective mass wN however is safe to freeze per substep
		m3d_vec3 rA = m3d_rotate(a.rotation, mp.localAnchorA);
		m3d_vec3 rB = m3d_rotate(b.rotation, mp.localAnchorB);

		float sep = mp.separationOffset +
					m3d_dot(m3d_sub(m3d_add(b.position, rB), m3d_add(a.position, rA)), n);
		if (sep >= 0.0f)
		{
			continue; // inactive (speculative point not yet reached)
		}

		// cap the depenetration per substep: an uncapped projection acts like
		// an impulsive kick when a rocking edge digs in, pumping energy
		if (sep < -maxCorrection)
		{
			sep = -maxCorrection;
		}

		float dLambda = -sep / mp.wN;
		mp.lambdaN += dLambda;
		ApplyAxisCorrection(b, dLambda, n, mp.angJB);
		ApplyAxisCorrection(a, -dLambda, n, mp.angJA);
	}
}

// Static friction runs ONCE per substep after all normal iterations. The
// anchors of a cached manifold are persistent stiction targets: at creation
// the anchor pair is tangentially coincident, so the tangential part of the
// current anchor gap IS the accumulated drift. Correcting toward this fixed
// target converges geometrically even under-relaxed - no creep accumulation
// (a per-substep re-anchored target integrates drift instead and either
// creeps or, at full gain, pumps tower sway).
static void ApplyContactFriction(m3d_world* world, Contact& c)
{
	Body& a = world->bodies[c.bodyA];
	Body& b = world->bodies[c.bodyB];
	Manifold& man = c.manifold;
	m3d_vec3 n = man.normal;

	for (int32_t i = 0; i < man.pointCount; ++i)
	{
		ManifoldPoint& mp = man.points[i];
		if (mp.lambdaN <= 0.0f)
		{
			continue;
		}
		m3d_vec3 rA = m3d_rotate(a.rotation, mp.frictionAnchorA);
		m3d_vec3 rB = m3d_rotate(b.rotation, mp.frictionAnchorB);
		m3d_vec3 gap = m3d_sub(m3d_add(b.position, rB), m3d_add(a.position, rA));
		m3d_vec3 dpt = m3d_sub(gap, m3d_scale(m3d_dot(gap, n), n));
		float lt = m3d_length(dpt);
		if (lt <= 1e-9f)
		{
			continue;
		}
		m3d_vec3 t = m3d_scale(1.0f / lt, dpt);
		float wt = GeneralizedInvMass(a, rA, t) + GeneralizedInvMass(b, rB, t);
		if (wt <= 0.0f)
		{
			continue;
		}
		float needed = kFrictionRelaxation * (lt / wt);
		float allowed = c.friction * mp.lambdaN - mp.lambdaT;
		if (allowed <= 0.0f)
		{
			continue;
		}
		float dL = needed < allowed ? needed : allowed;
		mp.lambdaT += dL;
		m3d_vec3 pt = m3d_scale(dL, t);
		ApplyPositionCorrection(b, m3d_neg(pt), rB);
		ApplyPositionCorrection(a, pt, rA);

		// the cone cannot hold the stretch: the pair is sliding, so the
		// stiction targets re-anchor at the current contact points
		if (needed > 2.0f * allowed)
		{
			mp.frictionAnchorA = mp.localAnchorA;
			mp.frictionAnchorB = mp.localAnchorB;
		}
	}
}

// ---------------------------------------------------------------------------
// Contact velocity pass: dynamic friction + restitution
// ---------------------------------------------------------------------------

static void SolveContactVelocity(m3d_world* world, Contact& c, float h, float gravityLen)
{
	Body& a = world->bodies[c.bodyA];
	Body& b = world->bodies[c.bodyB];
	Manifold& man = c.manifold;
	m3d_vec3 n = man.normal;
	float invH = h > 0.0f ? 1.0f / h : 0.0f;

	for (int32_t i = 0; i < man.pointCount; ++i)
	{
		ManifoldPoint& mp = man.points[i];
		if (mp.lambdaN <= 0.0f)
		{
			continue;
		}

		m3d_vec3 vrel = m3d_sub(m3d_add(b.linearVelocity, m3d_cross(b.angularVelocity, mp.rB)),
								m3d_add(a.linearVelocity, m3d_cross(a.angularVelocity, mp.rA)));
		float vn = m3d_dot(vrel, n);
		m3d_vec3 vt = m3d_sub(vrel, m3d_scale(vn, n));
		float vtLen = m3d_length(vt);

		// dynamic friction: decelerate tangential motion, capped by the
		// Coulomb limit derived from the positional normal correction
		if (vtLen > 1e-9f)
		{
			m3d_vec3 t = m3d_scale(1.0f / vtLen, vt);
			float wt = GeneralizedInvMass(a, mp.rA, t) + GeneralizedInvMass(b, mp.rB, t);
			if (wt > 0.0f)
			{
				float dv = c.friction * mp.lambdaN * invH; // = h * mu * |fn|, fn = lambdaN / h^2
				dv = dv < vtLen ? dv : vtLen;
				m3d_vec3 p = m3d_scale(-dv / wt, t);
				ApplyVelocityImpulse(b, p, mp.rB);
				ApplyVelocityImpulse(a, m3d_neg(p), mp.rA);
			}
		}

		// restitution against the pre-substep normal velocity
		m3d_vec3 vrelPrev = m3d_sub(m3d_add(b.prevLinearVelocity, m3d_cross(b.prevAngularVelocity, mp.rB)),
									m3d_add(a.prevLinearVelocity, m3d_cross(a.prevAngularVelocity, mp.rA)));
		float vnPrev = m3d_dot(vrelPrev, n);
		float e = fabsf(vnPrev) > 2.0f * gravityLen * h ? c.restitution : 0.0f;

		vrel = m3d_sub(m3d_add(b.linearVelocity, m3d_cross(b.angularVelocity, mp.rB)),
					   m3d_add(a.linearVelocity, m3d_cross(a.angularVelocity, mp.rA)));
		vn = m3d_dot(vrel, n);
		float target = -e * vnPrev;
		target = target > 0.0f ? target : 0.0f;
		// applied unconditionally (paper eq. 29): this also absorbs the
		// fictitious separating velocity created by position projection,
		// otherwise deep stacks chatter and never fall asleep
		if (mp.wN > 0.0f)
		{
			float j = (target - vn) / mp.wN;
			ApplyAxisVelocity(b, j, n, mp.angJB);
			ApplyAxisVelocity(a, -j, n, mp.angJA);
		}
	}
}

// ---------------------------------------------------------------------------
// Joints
// ---------------------------------------------------------------------------

static void PrepareJoint(m3d_world* world, Joint& j)
{
	j.compliance = 0.0f;
	if (j.type == JointType::Distance && j.hertz > 0.0f)
	{
		Body& a = world->bodies[j.bodyA];
		Body& b = world->bodies[j.bodyB];
		m3d_vec3 rA = m3d_rotate(a.rotation, j.localAnchorA);
		m3d_vec3 rB = m3d_rotate(b.rotation, j.localAnchorB);
		m3d_vec3 d = m3d_sub(m3d_add(b.position, rB), m3d_add(a.position, rA));
		m3d_vec3 u = m3d_normalize(d);
		float w = GeneralizedInvMass(a, rA, u) + GeneralizedInvMass(b, rB, u);
		if (w > 0.0f)
		{
			float meff = 1.0f / w;
			float omega = 2.0f * M3D_PI * j.hertz;
			float k = meff * omega * omega;
			j.compliance = k > 0.0f ? 1.0f / k : 0.0f;
		}
	}
}

static void SolveJointPosition(m3d_world* world, Joint& j, float h)
{
	Body& a = world->bodies[j.bodyA];
	Body& b = world->bodies[j.bodyB];
	m3d_vec3 rA = m3d_rotate(a.rotation, j.localAnchorA);
	m3d_vec3 rB = m3d_rotate(b.rotation, j.localAnchorB);
	m3d_vec3 d = m3d_sub(m3d_add(b.position, rB), m3d_add(a.position, rA));

	float len = m3d_length(d);
	if (len < 1e-9f)
	{
		return;
	}
	m3d_vec3 u = m3d_scale(1.0f / len, d);

	float C;
	if (j.type == JointType::Distance)
	{
		C = len - j.length;
	}
	else // Ball: drive the anchor gap to zero along its own direction
	{
		C = len;
	}

	float w = GeneralizedInvMass(a, rA, u) + GeneralizedInvMass(b, rB, u);
	float alphaTilde = j.compliance / (h * h);
	float denom = w + alphaTilde;
	if (denom <= 0.0f)
	{
		return;
	}
	float dLambda = (-C - alphaTilde * j.lambda) / denom;
	j.lambda += dLambda;

	m3d_vec3 p = m3d_scale(dLambda, u);
	ApplyPositionCorrection(b, p, rB);
	ApplyPositionCorrection(a, m3d_neg(p), rA);
}

// ---------------------------------------------------------------------------
// Per-substep stage helpers (shared by the serial and the graph-colored path)
// ---------------------------------------------------------------------------

static void IntegrateBodySubstep(m3d_world* world, int32_t bi, float h, m3d_vec3 gravity)
{
	Body& body = world->bodies[bi];
	if (body.type != M3D_BODY_DYNAMIC)
	{
		return;
	}
	body.prevLinearVelocity = body.linearVelocity;
	body.prevAngularVelocity = body.angularVelocity;

	body.linearVelocity = m3d_mul_add(body.linearVelocity, h * body.gravityScale, gravity);
	body.linearVelocity = m3d_scale(1.0f / (1.0f + h * body.linearDamping), body.linearVelocity);
	body.angularVelocity = m3d_scale(1.0f / (1.0f + h * body.angularDamping), body.angularVelocity);

	body.prevPosition = body.position;
	body.prevRotation = body.rotation;
	body.position = m3d_mul_add(body.position, h, body.linearVelocity);
	body.rotation = m3d_quat_integrate(body.rotation, body.angularVelocity, h);
	UpdateInvInertiaWorld(body);
}

// narrowphase at current substep poses, cached while bodies are still:
// resting contacts skip collision entirely and keep their anchors as
// persistent stiction targets for static friction
static void PrepareContactSubstep(m3d_world* world, int32_t ci, float h)
{
	Contact& c = world->contacts[ci];
	Body& a = world->bodies[c.bodyA];
	Body& b = world->bodies[c.bodyB];
	bool regen = !c.hasCache || PoseChanged(a, c.cachePosA, c.cacheRotA) || PoseChanged(b, c.cachePosB, c.cacheRotB);
	if (regen)
	{
		Manifold oldManifold = c.manifold;
		float relSpeed = m3d_length(m3d_sub(a.linearVelocity, b.linearVelocity));
		CollideShapes(a.shape, BodyTransform(a), b.shape, BodyTransform(b), kSpeculativeBase + h * relSpeed,
					  &c.manifold);
		InheritFrictionAnchors(oldManifold, c.manifold, a, b);
		c.hasCache = true;
		c.cachePosA = a.position;
		c.cacheRotA = a.rotation;
		c.cachePosB = b.position;
		c.cacheRotB = b.rotation;
	}
	// freeze anchors, the normal mass and the angular Jacobians for this substep
	m3d_vec3 n = c.manifold.normal;
	for (int32_t i = 0; i < c.manifold.pointCount; ++i)
	{
		ManifoldPoint& mp = c.manifold.points[i];
		mp.lambdaN = 0.0f;
		mp.lambdaT = 0.0f;
		mp.rA = m3d_rotate(a.rotation, mp.localAnchorA);
		mp.rB = m3d_rotate(b.rotation, mp.localAnchorB);
		m3d_vec3 crossA = m3d_cross(mp.rA, n);
		m3d_vec3 crossB = m3d_cross(mp.rB, n);
		mp.angJA = m3d_mat3_mulv(&a.invInertiaWorld, crossA);
		mp.angJB = m3d_mat3_mulv(&b.invInertiaWorld, crossB);
		mp.wN = a.invMass + b.invMass + m3d_dot(crossA, mp.angJA) + m3d_dot(crossB, mp.angJB);
		if (regen)
		{
			float ref = m3d_dot(m3d_sub(m3d_add(b.position, mp.rB), m3d_add(a.position, mp.rA)), n);
			mp.separationOffset = mp.separation - ref;
		}
	}
}

static void ReconstructBodyVelocity(m3d_world* world, int32_t bi, float invH)
{
	Body& body = world->bodies[bi];
	if (body.type != M3D_BODY_DYNAMIC)
	{
		return;
	}
	body.linearVelocity = m3d_scale(invH, m3d_sub(body.position, body.prevPosition));
	m3d_quat dq = m3d_quat_mul(body.rotation, m3d_quat_conj(body.prevRotation));
	float sign = dq.w >= 0.0f ? 1.0f : -1.0f;
	body.angularVelocity = m3d_scale(2.0f * invH * sign, m3d_v3(dq.x, dq.y, dq.z));
}

// ---------------------------------------------------------------------------
// Island solve (serial Gauss-Seidel; islands run in parallel with each other)
// ---------------------------------------------------------------------------

void SolveIslandXPBD(m3d_world* world, const std::vector<int32_t>& bodyIndices,
					 const std::vector<int32_t>& contactIndices, const std::vector<int32_t>& jointIndices, float dt,
					 int substepCount)
{
	float h = dt / (float)substepCount;
	float invH = h > 0.0f ? 1.0f / h : 0.0f;
	m3d_vec3 gravity = world->def.gravity;
	float gravityLen = m3d_length(gravity);
	int32_t positionIterations = world->def.positionIterations > 0 ? world->def.positionIterations : 2;

	for (int32_t bi : bodyIndices)
	{
		UpdateInvInertiaWorld(world->bodies[bi]);
	}
	for (int32_t ji : jointIndices)
	{
		PrepareJoint(world, world->joints[ji]);
	}

	for (int sub = 0; sub < substepCount; ++sub)
	{
		for (int32_t bi : bodyIndices)
		{
			IntegrateBodySubstep(world, bi, h, gravity);
		}
		for (int32_t ci : contactIndices)
		{
			PrepareContactSubstep(world, ci, h);
		}
		for (int32_t ji : jointIndices)
		{
			world->joints[ji].lambda = 0.0f;
		}

		// position projection; depenetration capped at kMaxDepenetrationSpeed
		float maxCorrection = kMaxDepenetrationSpeed * h;
		for (int32_t iter = 0; iter < positionIterations; ++iter)
		{
			for (int32_t ji : jointIndices)
			{
				SolveJointPosition(world, world->joints[ji], h);
			}
			for (int32_t ci : contactIndices)
			{
				SolveContactNormals(world, world->contacts[ci], maxCorrection);
			}
		}

		// friction once per substep, after the normals have converged
		for (int32_t ci : contactIndices)
		{
			ApplyContactFriction(world, world->contacts[ci]);
		}

		for (int32_t bi : bodyIndices)
		{
			ReconstructBodyVelocity(world, bi, invH);
		}

		// velocity pass
		for (int32_t ci : contactIndices)
		{
			SolveContactVelocity(world, world->contacts[ci], h, gravityLen);
		}
	}
}

// ---------------------------------------------------------------------------
// SIMD contact packets for the colored normal projection
// ---------------------------------------------------------------------------

static void BuildColorPackets(m3d_world* world, const std::vector<int32_t>& bucket,
							  std::vector<ContactPacket8>& packs)
{
	int32_t count = (int32_t)bucket.size();
	int32_t packetCount = (count + 7) / 8;
	packs.resize(packetCount);

	for (int32_t p = 0; p < packetCount; ++p)
	{
		ContactPacket8& pk = packs[p];
		pk.maxPoints = 0;
		for (int32_t l = 0; l < 8; ++l)
		{
			int32_t k = p * 8 + l;
			bool pad = k >= count;
			int32_t ci = bucket[pad ? p * 8 : k]; // padding copies lane 0 (safe loads)
			const Contact& c = world->contacts[ci];
			const Body& a = world->bodies[c.bodyA];
			const Body& b = world->bodies[c.bodyB];

			pk.contactIndex[l] = pad ? -1 : ci;
			pk.bodyA[l] = c.bodyA;
			pk.bodyB[l] = c.bodyB;
			pk.invMassA[l] = a.invMass;
			pk.invMassB[l] = b.invMass;
			pk.nx[l] = c.manifold.normal.x;
			pk.ny[l] = c.manifold.normal.y;
			pk.nz[l] = c.manifold.normal.z;

			int32_t pts = pad ? 0 : c.manifold.pointCount;
			pk.lanePoints[l] = pts;
			pk.maxPoints = pts > pk.maxPoints ? pts : pk.maxPoints;
			for (int32_t j = 0; j < 4; ++j)
			{
				bool on = j < pts;
				const ManifoldPoint& mp = c.manifold.points[on ? j : 0];
				pk.active[j][l] = on ? 1.0f : 0.0f;
				pk.lax[j][l] = mp.localAnchorA.x;
				pk.lay[j][l] = mp.localAnchorA.y;
				pk.laz[j][l] = mp.localAnchorA.z;
				pk.lbx[j][l] = mp.localAnchorB.x;
				pk.lby[j][l] = mp.localAnchorB.y;
				pk.lbz[j][l] = mp.localAnchorB.z;
				pk.sepOff[j][l] = mp.separationOffset;
				pk.invW[j][l] = (on && mp.wN > 0.0f) ? 1.0f / mp.wN : 1.0f;
				pk.jAx[j][l] = mp.angJA.x;
				pk.jAy[j][l] = mp.angJA.y;
				pk.jAz[j][l] = mp.angJA.z;
				pk.jBx[j][l] = mp.angJB.x;
				pk.jBy[j][l] = mp.angJB.y;
				pk.jBz[j][l] = mp.angJB.z;
				pk.lamN[j][l] = 0.0f;
			}
		}
	}
}

// One Gauss-Seidel sweep over a packet: 8 body-disjoint contacts in lanes,
// manifold points in sequence. Body state lives in lane arrays for the
// whole sweep; every inner loop is branchless and 8-wide.
static void SolvePacketNormals(m3d_world* world, ContactPacket8& pk, float maxCorrection)
{
	alignas(32) float pax[8], pay[8], paz[8], qax[8], qay[8], qaz[8], qaw[8];
	alignas(32) float pbx[8], pby[8], pbz[8], qbx[8], qby[8], qbz[8], qbw[8];

	for (int32_t l = 0; l < 8; ++l)
	{
		const Body& a = world->bodies[pk.bodyA[l]];
		const Body& b = world->bodies[pk.bodyB[l]];
		pax[l] = a.position.x; pay[l] = a.position.y; paz[l] = a.position.z;
		qax[l] = a.rotation.x; qay[l] = a.rotation.y; qaz[l] = a.rotation.z; qaw[l] = a.rotation.w;
		pbx[l] = b.position.x; pby[l] = b.position.y; pbz[l] = b.position.z;
		qbx[l] = b.rotation.x; qby[l] = b.rotation.y; qbz[l] = b.rotation.z; qbw[l] = b.rotation.w;
	}

	for (int32_t j = 0; j < pk.maxPoints; ++j)
	{
		for (int32_t l = 0; l < 8; ++l)
		{
			// rotate local anchors by the current rotations (r = v + w*t + q x t, t = 2 q x v)
			float vx = pk.lax[j][l], vy = pk.lay[j][l], vz = pk.laz[j][l];
			float tx = 2.0f * (qay[l] * vz - qaz[l] * vy);
			float ty = 2.0f * (qaz[l] * vx - qax[l] * vz);
			float tz = 2.0f * (qax[l] * vy - qay[l] * vx);
			float rax = vx + qaw[l] * tx + (qay[l] * tz - qaz[l] * ty);
			float ray = vy + qaw[l] * ty + (qaz[l] * tx - qax[l] * tz);
			float raz = vz + qaw[l] * tz + (qax[l] * ty - qay[l] * tx);

			vx = pk.lbx[j][l]; vy = pk.lby[j][l]; vz = pk.lbz[j][l];
			tx = 2.0f * (qby[l] * vz - qbz[l] * vy);
			ty = 2.0f * (qbz[l] * vx - qbx[l] * vz);
			tz = 2.0f * (qbx[l] * vy - qby[l] * vx);
			float rbx = vx + qbw[l] * tx + (qby[l] * tz - qbz[l] * ty);
			float rby = vy + qbw[l] * ty + (qbz[l] * tx - qbx[l] * tz);
			float rbz = vz + qbw[l] * tz + (qbx[l] * ty - qby[l] * tx);

			float sep = pk.sepOff[j][l] + (pbx[l] + rbx - pax[l] - rax) * pk.nx[l] +
						(pby[l] + rby - pay[l] - ray) * pk.ny[l] + (pbz[l] + rbz - paz[l] - raz) * pk.nz[l];

			// branchless: active points with sep < 0 only, capped depenetration
			float gate = sep < 0.0f ? pk.active[j][l] : 0.0f;
			float capped = sep < -maxCorrection ? -maxCorrection : sep;
			float dL = -capped * pk.invW[j][l] * gate;
			pk.lamN[j][l] += dL;

			// linear response
			pbx[l] += dL * pk.invMassB[l] * pk.nx[l];
			pby[l] += dL * pk.invMassB[l] * pk.ny[l];
			pbz[l] += dL * pk.invMassB[l] * pk.nz[l];
			pax[l] -= dL * pk.invMassA[l] * pk.nx[l];
			pay[l] -= dL * pk.invMassA[l] * pk.ny[l];
			paz[l] -= dL * pk.invMassA[l] * pk.nz[l];

			// angular response via frozen Jacobians: q += 0.5*(dw x q), fast renorm
			float dwx = dL * pk.jBx[j][l], dwy = dL * pk.jBy[j][l], dwz = dL * pk.jBz[j][l];
			float dqx = dwx * qbw[l] + dwy * qbz[l] - dwz * qby[l];
			float dqy = -dwx * qbz[l] + dwy * qbw[l] + dwz * qbx[l];
			float dqz = dwx * qby[l] - dwy * qbx[l] + dwz * qbw[l];
			float dqw = -(dwx * qbx[l] + dwy * qby[l] + dwz * qbz[l]);
			float x = qbx[l] + 0.5f * dqx, y = qby[l] + 0.5f * dqy, z = qbz[l] + 0.5f * dqz, w = qbw[l] + 0.5f * dqw;
			float s = 1.5f - 0.5f * (x * x + y * y + z * z + w * w);
			qbx[l] = x * s; qby[l] = y * s; qbz[l] = z * s; qbw[l] = w * s;

			dwx = -dL * pk.jAx[j][l]; dwy = -dL * pk.jAy[j][l]; dwz = -dL * pk.jAz[j][l];
			dqx = dwx * qaw[l] + dwy * qaz[l] - dwz * qay[l];
			dqy = -dwx * qaz[l] + dwy * qaw[l] + dwz * qax[l];
			dqz = dwx * qay[l] - dwy * qax[l] + dwz * qaw[l];
			dqw = -(dwx * qax[l] + dwy * qay[l] + dwz * qaz[l]);
			x = qax[l] + 0.5f * dqx; y = qay[l] + 0.5f * dqy; z = qaz[l] + 0.5f * dqz; w = qaw[l] + 0.5f * dqw;
			s = 1.5f - 0.5f * (x * x + y * y + z * z + w * w);
			qax[l] = x * s; qay[l] = y * s; qaz[l] = z * s; qaw[l] = w * s;
		}
	}

	// write back dynamic bodies of real lanes (statics are lane-shared, read-only)
	for (int32_t l = 0; l < 8; ++l)
	{
		if (pk.contactIndex[l] < 0)
		{
			continue;
		}
		if (pk.invMassA[l] > 0.0f)
		{
			Body& a = world->bodies[pk.bodyA[l]];
			a.position = m3d_v3(pax[l], pay[l], paz[l]);
			a.rotation = { qax[l], qay[l], qaz[l], qaw[l] };
		}
		if (pk.invMassB[l] > 0.0f)
		{
			Body& b = world->bodies[pk.bodyB[l]];
			b.position = m3d_v3(pbx[l], pby[l], pbz[l]);
			b.rotation = { qbx[l], qby[l], qbz[l], qbw[l] };
		}
	}
}

static void FlushPacketLambdas(m3d_world* world, const std::vector<ContactPacket8>& packs)
{
	for (const ContactPacket8& pk : packs)
	{
		for (int32_t l = 0; l < 8; ++l)
		{
			int32_t ci = pk.contactIndex[l];
			if (ci < 0)
			{
				continue;
			}
			Contact& c = world->contacts[ci];
			for (int32_t j = 0; j < pk.lanePoints[l]; ++j)
			{
				c.manifold.points[j].lambdaN = pk.lamN[j][l];
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Graph-colored island solve: parallelism INSIDE one big island
// ---------------------------------------------------------------------------

void SolveIslandColoredXPBD(m3d_world* world, const std::vector<int32_t>& bodyIndices,
							const std::vector<int32_t>& contactIndices, const std::vector<int32_t>& jointIndices,
							float dt, int substepCount)
{
	float h = dt / (float)substepCount;
	float invH = h > 0.0f ? 1.0f / h : 0.0f;
	m3d_vec3 gravity = world->def.gravity;
	float gravityLen = m3d_length(gravity);
	int32_t positionIterations = world->def.positionIterations > 0 ? world->def.positionIterations : 2;
	ThreadPool* pool = world->pool;
	int32_t workers = pool->WorkerCount();

	// --- greedy coloring: no two contacts in a bucket share a dynamic body.
	// Static/kinematic bodies are read-only in the solve and take no color,
	// so hundreds of ground contacts spread across all buckets.
	if (world->colorMask.size() < world->bodies.size())
	{
		world->colorMask.resize(world->bodies.size(), 0);
	}
	for (int32_t bi : bodyIndices)
	{
		world->colorMask[bi] = 0;
	}
	if ((int32_t)world->colorBuckets.size() < kMaxGraphColors)
	{
		world->colorBuckets.resize(kMaxGraphColors);
	}
	int32_t colorCount = 0;
	for (int32_t c = 0; c < kMaxGraphColors; ++c)
	{
		world->colorBuckets[c].clear();
	}
	world->overflowContacts.clear();

	for (int32_t ci : contactIndices)
	{
		Contact& c = world->contacts[ci];
		Body& a = world->bodies[c.bodyA];
		Body& b = world->bodies[c.bodyB];
		bool dynA = a.invMass > 0.0f;
		bool dynB = b.invMass > 0.0f;
		uint32_t used = 0;
		if (dynA)
		{
			used |= world->colorMask[c.bodyA];
		}
		if (dynB)
		{
			used |= world->colorMask[c.bodyB];
		}
		uint32_t freeBits = ~used;
		int32_t color = freeBits != 0 ? __builtin_ctz(freeBits) : kMaxGraphColors;
		if (color >= kMaxGraphColors)
		{
			world->overflowContacts.push_back(ci);
			continue;
		}
		world->colorBuckets[color].push_back(ci);
		uint32_t bit = 1u << color;
		if (dynA)
		{
			world->colorMask[c.bodyA] |= bit;
		}
		if (dynB)
		{
			world->colorMask[c.bodyB] |= bit;
		}
		colorCount = color + 1 > colorCount ? color + 1 : colorCount;
	}

	auto colorGrain = [&](int32_t size) {
		int32_t g = size / (3 * workers) + 1;
		return g < 4 ? 4 : g;
	};

	const int32_t bodyCount = (int32_t)bodyIndices.size();
	const int32_t contactCount = (int32_t)contactIndices.size();

	pool->ParallelFor(bodyCount, 128, [&](int32_t k) { UpdateInvInertiaWorld(world->bodies[bodyIndices[k]]); });
	for (int32_t ji : jointIndices)
	{
		PrepareJoint(world, world->joints[ji]);
	}

	if ((int32_t)world->colorPackets.size() < kMaxGraphColors)
	{
		world->colorPackets.resize(kMaxGraphColors);
	}

	for (int sub = 0; sub < substepCount; ++sub)
	{
		pool->ParallelFor(bodyCount, 64, [&](int32_t k) { IntegrateBodySubstep(world, bodyIndices[k], h, gravity); });
		pool->ParallelFor(contactCount, 32, [&](int32_t k) { PrepareContactSubstep(world, contactIndices[k], h); });
		for (int32_t ji : jointIndices)
		{
			world->joints[ji].lambda = 0.0f;
		}

		// pack the color buckets into 8-wide SIMD packets for this substep
		for (int32_t color = 0; color < colorCount; ++color)
		{
			BuildColorPackets(world, world->colorBuckets[color], world->colorPackets[color]);
		}

		float maxCorrection = kMaxDepenetrationSpeed * h;
		for (int32_t iter = 0; iter < positionIterations; ++iter)
		{
			for (int32_t ji : jointIndices)
			{
				SolveJointPosition(world, world->joints[ji], h);
			}
			for (int32_t color = 0; color < colorCount; ++color)
			{
				auto& packs = world->colorPackets[color];
				pool->ParallelFor((int32_t)packs.size(), 1,
								  [&](int32_t k) { SolvePacketNormals(world, packs[k], maxCorrection); });
			}
			for (int32_t ci : world->overflowContacts)
			{
				SolveContactNormals(world, world->contacts[ci], maxCorrection);
			}
		}

		// hand the accumulated normal lambdas back to the manifold points:
		// friction cones and the velocity pass read them
		for (int32_t color = 0; color < colorCount; ++color)
		{
			FlushPacketLambdas(world, world->colorPackets[color]);
		}

		for (int32_t color = 0; color < colorCount; ++color)
		{
			auto& bucket = world->colorBuckets[color];
			pool->ParallelFor((int32_t)bucket.size(), colorGrain((int32_t)bucket.size()),
							  [&](int32_t k) { ApplyContactFriction(world, world->contacts[bucket[k]]); });
		}
		for (int32_t ci : world->overflowContacts)
		{
			ApplyContactFriction(world, world->contacts[ci]);
		}

		pool->ParallelFor(bodyCount, 64, [&](int32_t k) { ReconstructBodyVelocity(world, bodyIndices[k], invH); });

		for (int32_t color = 0; color < colorCount; ++color)
		{
			auto& bucket = world->colorBuckets[color];
			pool->ParallelFor((int32_t)bucket.size(), colorGrain((int32_t)bucket.size()),
							  [&](int32_t k) { SolveContactVelocity(world, world->contacts[bucket[k]], h, gravityLen); });
		}
		for (int32_t ci : world->overflowContacts)
		{
			SolveContactVelocity(world, world->contacts[ci], h, gravityLen);
		}
	}
}

} // namespace m3d
