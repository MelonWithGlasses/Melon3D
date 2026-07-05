// Narrowphase collision: contact manifold generation and shape ray casts.
// Manifold normals point from shape A to shape B. All contacts may be
// speculative: points with positive separation up to margin
// are kept so the solver can prevent tunneling.
#include "m3d_internal.h"

#include <float.h>

namespace m3d
{

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static void AddManifoldPoint(Manifold* m, m3d_transform xfA, m3d_transform xfB, m3d_vec3 pointOnA, m3d_vec3 pointOnB,
							 float separation, uint32_t id)
{
	if (m->pointCount >= 4)
	{
		return;
	}
	ManifoldPoint& mp = m->points[m->pointCount++];
	mp.localAnchorA = m3d_inv_transform_point(xfA, pointOnA);
	mp.localAnchorB = m3d_inv_transform_point(xfB, pointOnB);
	mp.frictionAnchorA = mp.localAnchorA;
	mp.frictionAnchorB = mp.localAnchorB;
	mp.separation = separation;
	mp.id = id;
	mp.lambdaN = 0.0f;
	mp.lambdaT = 0.0f;
	mp.separationOffset = 0.0f;
}

static m3d_vec3 ClosestPointOnSegment(m3d_vec3 p, m3d_vec3 a, m3d_vec3 b)
{
	m3d_vec3 ab = m3d_sub(b, a);
	float denom = m3d_dot(ab, ab);
	if (denom < 1e-12f)
	{
		return a;
	}
	float t = m3d_dot(m3d_sub(p, a), ab) / denom;
	t = t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);
	return m3d_mul_add(a, t, ab);
}

// Ericson, Real-Time Collision Detection 5.1.9
static void ClosestPtSegmentSegment(m3d_vec3 p1, m3d_vec3 q1, m3d_vec3 p2, m3d_vec3 q2, m3d_vec3* c1, m3d_vec3* c2)
{
	m3d_vec3 d1 = m3d_sub(q1, p1);
	m3d_vec3 d2 = m3d_sub(q2, p2);
	m3d_vec3 r = m3d_sub(p1, p2);
	float a = m3d_dot(d1, d1);
	float e = m3d_dot(d2, d2);
	float f = m3d_dot(d2, r);
	float s, t;

	if (a < 1e-12f && e < 1e-12f)
	{
		*c1 = p1;
		*c2 = p2;
		return;
	}
	if (a < 1e-12f)
	{
		s = 0.0f;
		t = f / e;
	}
	else
	{
		float c = m3d_dot(d1, r);
		if (e < 1e-12f)
		{
			t = 0.0f;
			s = -c / a;
		}
		else
		{
			float b = m3d_dot(d1, d2);
			float denom = a * e - b * b;
			s = denom > 1e-12f ? (b * f - c * e) / denom : 0.0f;
			s = s < 0.0f ? 0.0f : (s > 1.0f ? 1.0f : s);
			t = (b * s + f) / e;
			if (t < 0.0f)
			{
				t = 0.0f;
				s = -c / a;
			}
			else if (t > 1.0f)
			{
				t = 1.0f;
				s = (b - c) / a;
			}
			s = s < 0.0f ? 0.0f : (s > 1.0f ? 1.0f : s);
		}
	}
	*c1 = m3d_mul_add(p1, s, d1);
	*c2 = m3d_mul_add(p2, t, d2);
}

static m3d_vec3 ClampToBox(m3d_vec3 p, m3d_vec3 e)
{
	m3d_vec3 q;
	q.x = p.x < -e.x ? -e.x : (p.x > e.x ? e.x : p.x);
	q.y = p.y < -e.y ? -e.y : (p.y > e.y ? e.y : p.y);
	q.z = p.z < -e.z ? -e.z : (p.z > e.z ? e.z : p.z);
	return q;
}

// ---------------------------------------------------------------------------
// Sphere - Sphere
// ---------------------------------------------------------------------------

static void CollideSpheres(const Shape& a, m3d_transform xfA, const Shape& b, m3d_transform xfB, float margin, Manifold* m)
{
	m3d_vec3 d = m3d_sub(xfB.p, xfA.p);
	float dist = m3d_length(d);
	m3d_vec3 n = dist > 1e-9f ? m3d_scale(1.0f / dist, d) : m3d_v3(0.0f, 1.0f, 0.0f);
	float separation = dist - a.radius - b.radius;
	if (separation > margin)
	{
		return;
	}
	m->normal = n;
	m3d_vec3 pA = m3d_mul_add(xfA.p, a.radius, n);
	m3d_vec3 pB = m3d_mul_add(xfB.p, -b.radius, n);
	AddManifoldPoint(m, xfA, xfB, pA, pB, separation, 0);
}

// ---------------------------------------------------------------------------
// Sphere (A) - Capsule (B)
// ---------------------------------------------------------------------------

static void CollideSphereCapsule(const Shape& a, m3d_transform xfA, const Shape& b, m3d_transform xfB, float margin, Manifold* m)
{
	m3d_vec3 axis = m3d_rotate(xfB.q, m3d_v3(0.0f, b.halfHeight, 0.0f));
	m3d_vec3 p1 = m3d_sub(xfB.p, axis);
	m3d_vec3 p2 = m3d_add(xfB.p, axis);
	m3d_vec3 c = ClosestPointOnSegment(xfA.p, p1, p2);

	m3d_vec3 d = m3d_sub(c, xfA.p);
	float dist = m3d_length(d);
	m3d_vec3 n = dist > 1e-9f ? m3d_scale(1.0f / dist, d) : m3d_v3(0.0f, 1.0f, 0.0f);
	float separation = dist - a.radius - b.radius;
	if (separation > margin)
	{
		return;
	}
	m->normal = n;
	m3d_vec3 pA = m3d_mul_add(xfA.p, a.radius, n);
	m3d_vec3 pB = m3d_mul_add(c, -b.radius, n);
	AddManifoldPoint(m, xfA, xfB, pA, pB, separation, 0);
}

// ---------------------------------------------------------------------------
// Capsule (A) - Capsule (B)
// ---------------------------------------------------------------------------

static void CollideCapsules(const Shape& a, m3d_transform xfA, const Shape& b, m3d_transform xfB, float margin, Manifold* m)
{
	m3d_vec3 axisA = m3d_rotate(xfA.q, m3d_v3(0.0f, a.halfHeight, 0.0f));
	m3d_vec3 a1 = m3d_sub(xfA.p, axisA);
	m3d_vec3 a2 = m3d_add(xfA.p, axisA);
	m3d_vec3 axisB = m3d_rotate(xfB.q, m3d_v3(0.0f, b.halfHeight, 0.0f));
	m3d_vec3 b1 = m3d_sub(xfB.p, axisB);
	m3d_vec3 b2 = m3d_add(xfB.p, axisB);

	m3d_vec3 uA = m3d_normalize(m3d_sub(a2, a1));
	m3d_vec3 uB = m3d_normalize(m3d_sub(b2, b1));
	float parallelism = m3d_length(m3d_cross(uA, uB));

	// Parallel capsules resting on each other need two contact points.
	if (parallelism < 0.05f && a.halfHeight > 1e-6f)
	{
		float lenA = 2.0f * a.halfHeight;
		float tb1 = m3d_dot(m3d_sub(b1, a1), uA);
		float tb2 = m3d_dot(m3d_sub(b2, a1), uA);
		float lo = fmaxf(0.0f, fminf(tb1, tb2));
		float hi = fminf(lenA, fmaxf(tb1, tb2));
		if (hi - lo > 1e-4f)
		{
			// shared normal from the interval midpoint
			m3d_vec3 midA = m3d_mul_add(a1, 0.5f * (lo + hi), uA);
			m3d_vec3 midB = ClosestPointOnSegment(midA, b1, b2);
			m3d_vec3 d = m3d_sub(midB, midA);
			float dist = m3d_length(d);
			m3d_vec3 n = dist > 1e-9f ? m3d_scale(1.0f / dist, d) : m3d_v3(0.0f, 1.0f, 0.0f);
			if (dist - a.radius - b.radius > margin)
			{
				return;
			}
			m->normal = n;
			float ts[2] = { lo, hi };
			for (int i = 0; i < 2; ++i)
			{
				m3d_vec3 cA = m3d_mul_add(a1, ts[i], uA);
				m3d_vec3 cB = ClosestPointOnSegment(cA, b1, b2);
				float di = m3d_length(m3d_sub(cB, cA));
				float sep = di - a.radius - b.radius;
				if (sep <= margin)
				{
					m3d_vec3 pA = m3d_mul_add(cA, a.radius, n);
					m3d_vec3 pB = m3d_mul_add(cB, -b.radius, n);
					AddManifoldPoint(m, xfA, xfB, pA, pB, sep, (uint32_t)i);
				}
			}
			if (m->pointCount > 0)
			{
				return;
			}
		}
	}

	m3d_vec3 cA, cB;
	ClosestPtSegmentSegment(a1, a2, b1, b2, &cA, &cB);
	m3d_vec3 d = m3d_sub(cB, cA);
	float dist = m3d_length(d);
	m3d_vec3 n = dist > 1e-9f ? m3d_scale(1.0f / dist, d) : m3d_v3(0.0f, 1.0f, 0.0f);
	float separation = dist - a.radius - b.radius;
	if (separation > margin)
	{
		return;
	}
	m->normal = n;
	m3d_vec3 pA = m3d_mul_add(cA, a.radius, n);
	m3d_vec3 pB = m3d_mul_add(cB, -b.radius, n);
	AddManifoldPoint(m, xfA, xfB, pA, pB, separation, 0);
}

// ---------------------------------------------------------------------------
// Sphere (A) - Box (B)
// ---------------------------------------------------------------------------

static void CollideSphereBox(const Shape& a, m3d_transform xfA, const Shape& b, m3d_transform xfB, float margin, Manifold* m)
{
	m3d_vec3 c = m3d_inv_transform_point(xfB, xfA.p);
	m3d_vec3 e = b.halfExtents;
	m3d_vec3 q = ClampToBox(c, e);
	m3d_vec3 d = m3d_sub(c, q);
	float dist2 = m3d_length_sq(d);

	m3d_vec3 nLocal; // from box toward sphere
	float separation;
	if (dist2 > 1e-12f)
	{
		float dist = sqrtf(dist2);
		nLocal = m3d_scale(1.0f / dist, d);
		separation = dist - a.radius;
	}
	else
	{
		// center inside the box: push out along the face of least penetration
		float dx = e.x - fabsf(c.x);
		float dy = e.y - fabsf(c.y);
		float dz = e.z - fabsf(c.z);
		if (dx < dy && dx < dz)
		{
			nLocal = m3d_v3(c.x >= 0.0f ? 1.0f : -1.0f, 0.0f, 0.0f);
			separation = -(dx + a.radius);
			q = m3d_v3(nLocal.x * e.x, c.y, c.z);
		}
		else if (dy < dz)
		{
			nLocal = m3d_v3(0.0f, c.y >= 0.0f ? 1.0f : -1.0f, 0.0f);
			separation = -(dy + a.radius);
			q = m3d_v3(c.x, nLocal.y * e.y, c.z);
		}
		else
		{
			nLocal = m3d_v3(0.0f, 0.0f, c.z >= 0.0f ? 1.0f : -1.0f);
			separation = -(dz + a.radius);
			q = m3d_v3(c.x, c.y, nLocal.z * e.z);
		}
	}

	if (separation > margin)
	{
		return;
	}

	m3d_vec3 nWorld = m3d_rotate(xfB.q, nLocal); // box -> sphere
	m->normal = m3d_neg(nWorld);                      // A (sphere) -> B (box)
	m3d_vec3 pB = m3d_transform_point(xfB, q);
	m3d_vec3 pA = m3d_mul_add(xfA.p, a.radius, m->normal);
	AddManifoldPoint(m, xfA, xfB, pA, pB, separation, 0);
}

// ---------------------------------------------------------------------------
// Capsule (A) - Box (B)
// ---------------------------------------------------------------------------

// Squared distance from a point to a box (local frame)
static float PointBoxDistSq(m3d_vec3 p, m3d_vec3 e)
{
	m3d_vec3 q = ClampToBox(p, e);
	return m3d_length_sq(m3d_sub(p, q));
}

static void CollideCapsuleBox(const Shape& a, m3d_transform xfA, const Shape& b, m3d_transform xfB, float margin, Manifold* m)
{
	m3d_vec3 e = b.halfExtents;
	float r = a.radius;

	m3d_vec3 axis = m3d_rotate(xfA.q, m3d_v3(0.0f, a.halfHeight, 0.0f));
	m3d_vec3 segA = m3d_inv_transform_point(xfB, m3d_sub(xfA.p, axis));
	m3d_vec3 segB = m3d_inv_transform_point(xfB, m3d_add(xfA.p, axis));
	m3d_vec3 segDir = m3d_sub(segB, segA);

	// Distance from a convex set along the segment is convex in t: ternary search.
	float lo = 0.0f, hi = 1.0f;
	for (int i = 0; i < 40; ++i)
	{
		float t1 = lo + (hi - lo) / 3.0f;
		float t2 = hi - (hi - lo) / 3.0f;
		float d1 = PointBoxDistSq(m3d_mul_add(segA, t1, segDir), e);
		float d2 = PointBoxDistSq(m3d_mul_add(segA, t2, segDir), e);
		if (d1 < d2)
		{
			hi = t2;
		}
		else
		{
			lo = t1;
		}
	}
	float tStar = 0.5f * (lo + hi);
	m3d_vec3 c = m3d_mul_add(segA, tStar, segDir);
	m3d_vec3 q = ClampToBox(c, e);
	m3d_vec3 diff = m3d_sub(c, q);
	float dist2 = m3d_length_sq(diff);

	m3d_vec3 nLocal; // box -> capsule
	float separation;
	bool deep = dist2 <= 1e-10f;
	if (!deep)
	{
		float dist = sqrtf(dist2);
		nLocal = m3d_scale(1.0f / dist, diff);
		separation = dist - r;
		if (separation > margin)
		{
			return;
		}
	}
	else
	{
		// Segment touches/penetrates the box: pick the face of least penetration
		float dx = e.x - fabsf(c.x);
		float dy = e.y - fabsf(c.y);
		float dz = e.z - fabsf(c.z);
		if (dx < dy && dx < dz)
		{
			nLocal = m3d_v3(c.x >= 0.0f ? 1.0f : -1.0f, 0.0f, 0.0f);
			separation = -(dx + r);
		}
		else if (dy < dz)
		{
			nLocal = m3d_v3(0.0f, c.y >= 0.0f ? 1.0f : -1.0f, 0.0f);
			separation = -(dy + r);
		}
		else
		{
			nLocal = m3d_v3(0.0f, 0.0f, c.z >= 0.0f ? 1.0f : -1.0f);
			separation = -(dz + r);
		}
		q = m3d_sub(c, m3d_scale(separation + r, nLocal)); // point on box surface
	}

	m3d_vec3 nWorldBoxToCap = m3d_rotate(xfB.q, nLocal);
	m->normal = m3d_neg(nWorldBoxToCap); // A (capsule) -> B (box)

	// Two contact points when the capsule lies parallel to the contact plane.
	m3d_vec3 u = m3d_normalize(segDir);
	if (fabsf(m3d_dot(u, nLocal)) < 0.05f && a.halfHeight > 1e-6f)
	{
		int added = 0;
		m3d_vec3 ends[2] = { segA, segB };
		for (int i = 0; i < 2; ++i)
		{
			m3d_vec3 p = ends[i];
			m3d_vec3 qp = ClampToBox(p, e);
			m3d_vec3 dp = m3d_sub(p, qp);
			float d2 = m3d_length_sq(dp);
			float sep;
			if (d2 > 1e-10f)
			{
				sep = sqrtf(d2) - r;
			}
			else
			{
				// endpoint inside the box: signed distance to the chosen face
				float extent = e.x * fabsf(nLocal.x) + e.y * fabsf(nLocal.y) + e.z * fabsf(nLocal.z);
				sep = (m3d_dot(p, nLocal) - extent) - r;
				qp = m3d_sub(p, m3d_scale(m3d_dot(p, nLocal) - extent, nLocal));
			}
			if (sep <= margin)
			{
				m3d_vec3 pCapWorld = m3d_mul_add(m3d_transform_point(xfB, p), r, m->normal);
				m3d_vec3 pBoxWorld = m3d_transform_point(xfB, qp);
				AddManifoldPoint(m, xfA, xfB, pCapWorld, pBoxWorld, sep, (uint32_t)i);
				++added;
			}
		}
		if (added > 0)
		{
			return;
		}
	}

	m3d_vec3 pCapWorld = m3d_mul_add(m3d_transform_point(xfB, c), r, m->normal);
	m3d_vec3 pBoxWorld = m3d_transform_point(xfB, q);
	AddManifoldPoint(m, xfA, xfB, pCapWorld, pBoxWorld, separation, 4);
}

// ---------------------------------------------------------------------------
// Box - Box (SAT with reference face clipping, edge-edge closest points)
// ---------------------------------------------------------------------------

struct ClipVertex
{
	m3d_vec3 v;
	uint32_t id;
};

// Sutherland-Hodgman clip of a polygon against the half-space dot(x, n) <= offset
static int ClipPolygon(const ClipVertex* in, int inCount, m3d_vec3 n, float offset, uint32_t planeId, ClipVertex* out)
{
	int outCount = 0;
	for (int i = 0; i < inCount; ++i)
	{
		const ClipVertex& v1 = in[i];
		const ClipVertex& v2 = in[(i + 1) % inCount];
		float d1 = m3d_dot(v1.v, n) - offset;
		float d2 = m3d_dot(v2.v, n) - offset;

		if (d1 <= 0.0f)
		{
			out[outCount++] = v1;
		}
		if (d1 * d2 < 0.0f && outCount < 8)
		{
			float t = d1 / (d1 - d2);
			out[outCount].v = m3d_lerp(v1.v, v2.v, t);
			out[outCount].id = ((v1.id ^ v2.id) << 8) | planeId;
			++outCount;
		}
		if (outCount >= 8)
		{
			break;
		}
	}
	return outCount;
}

// Support edge of an OBB in a given world direction
static void SupportEdge(m3d_transform xf, const m3d_mat3& rot, m3d_vec3 e, m3d_vec3 dirWorld, m3d_vec3* pOut, m3d_vec3* qOut)
{
	m3d_vec3 d = m3d_mat3_mulv_t(&rot, dirWorld); // to local
	m3d_vec3 absD = m3d_abs(d);

	int edgeAxis = 0;
	if (absD.y < absD.x && absD.y <= absD.z)
	{
		edgeAxis = 1;
	}
	else if (absD.z < absD.x && absD.z <= absD.y)
	{
		edgeAxis = 2;
	}
	else if (absD.x <= absD.y && absD.x <= absD.z)
	{
		edgeAxis = 0;
	}

	float sx = d.x >= 0.0f ? 1.0f : -1.0f;
	float sy = d.y >= 0.0f ? 1.0f : -1.0f;
	float sz = d.z >= 0.0f ? 1.0f : -1.0f;

	m3d_vec3 p = m3d_v3(sx * e.x, sy * e.y, sz * e.z);
	m3d_vec3 q = p;
	if (edgeAxis == 0)
	{
		p.x = -e.x;
		q.x = e.x;
	}
	else if (edgeAxis == 1)
	{
		p.y = -e.y;
		q.y = e.y;
	}
	else
	{
		p.z = -e.z;
		q.z = e.z;
	}
	*pOut = m3d_transform_point(xf, p);
	*qOut = m3d_transform_point(xf, q);
}

static void CollideBoxes(const Shape& a, m3d_transform xfA, const Shape& b, m3d_transform xfB, float margin, Manifold* m)
{
	m3d_mat3 rotA = m3d_quat_to_mat3(xfA.q);
	m3d_mat3 rotB = m3d_quat_to_mat3(xfB.q);
	m3d_vec3 eA = a.halfExtents;
	m3d_vec3 eB = b.halfExtents;
	m3d_vec3 d = m3d_sub(xfB.p, xfA.p);

	const m3d_vec3 axesA[3] = { rotA.cx, rotA.cy, rotA.cz };
	const m3d_vec3 axesB[3] = { rotB.cx, rotB.cy, rotB.cz };
	const float extA[3] = { eA.x, eA.y, eA.z };
	const float extB[3] = { eB.x, eB.y, eB.z };

	auto projectRadius = [](const m3d_vec3* axes, const float* ext, m3d_vec3 n) {
		return ext[0] * fabsf(m3d_dot(axes[0], n)) + ext[1] * fabsf(m3d_dot(axes[1], n)) + ext[2] * fabsf(m3d_dot(axes[2], n));
	};

	float faceASep = -FLT_MAX;
	int faceAIndex = 0;
	for (int i = 0; i < 3; ++i)
	{
		float s = fabsf(m3d_dot(d, axesA[i])) - (extA[i] + projectRadius(axesB, extB, axesA[i]));
		if (s > margin)
		{
			return;
		}
		if (s > faceASep)
		{
			faceASep = s;
			faceAIndex = i;
		}
	}

	float faceBSep = -FLT_MAX;
	int faceBIndex = 0;
	for (int i = 0; i < 3; ++i)
	{
		float s = fabsf(m3d_dot(d, axesB[i])) - (extB[i] + projectRadius(axesA, extA, axesB[i]));
		if (s > margin)
		{
			return;
		}
		if (s > faceBSep)
		{
			faceBSep = s;
			faceBIndex = i;
		}
	}

	float edgeSep = -FLT_MAX;
	int edgeI = 0, edgeJ = 0;
	m3d_vec3 edgeAxis = m3d_v3_zero();
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			m3d_vec3 axis = m3d_cross(axesA[i], axesB[j]);
			float len2 = m3d_length_sq(axis);
			if (len2 < 1e-8f)
			{
				continue; // near-parallel edges: face axes cover this direction
			}
			axis = m3d_scale(1.0f / sqrtf(len2), axis);
			float s = fabsf(m3d_dot(d, axis)) - (projectRadius(axesA, extA, axis) + projectRadius(axesB, extB, axis));
			if (s > margin)
			{
				return;
			}
			if (s > edgeSep)
			{
				edgeSep = s;
				edgeI = i;
				edgeJ = j;
				edgeAxis = axis;
			}
		}
	}

	// Prefer face contacts for coherence (qu3e-style tolerances)
	const float kRelTol = 0.95f;
	const float kAbsTol = 0.01f * kLinearSlop * 200.0f; // ~1cm scaled to slop
	float bestFaceSep = faceASep;
	bool refIsA = true;
	int refFace = faceAIndex;
	if (faceBSep > kRelTol * faceASep + kAbsTol)
	{
		bestFaceSep = faceBSep;
		refIsA = false;
		refFace = faceBIndex;
	}

	if (edgeSep > kRelTol * bestFaceSep + kAbsTol)
	{
		// Edge contact
		m3d_vec3 n = edgeAxis;
		if (m3d_dot(d, n) < 0.0f)
		{
			n = m3d_neg(n);
		}
		m3d_vec3 pA1, pA2, pB1, pB2;
		SupportEdge(xfA, rotA, eA, n, &pA1, &pA2);
		SupportEdge(xfB, rotB, eB, m3d_neg(n), &pB1, &pB2);
		m3d_vec3 cA, cB;
		ClosestPtSegmentSegment(pA1, pA2, pB1, pB2, &cA, &cB);
		m->normal = n;
		AddManifoldPoint(m, xfA, xfB, cA, cB, edgeSep, 0xE0000000u | (uint32_t)(edgeI * 3 + edgeJ));
		return;
	}

	// Face contact: clip the incident face of the other box against the
	// side planes of the reference face.
	m3d_transform xfRef = refIsA ? xfA : xfB;
	m3d_transform xfInc = refIsA ? xfB : xfA;
	const m3d_mat3& rotInc = refIsA ? rotB : rotA;
	const m3d_vec3* refAxes = refIsA ? axesA : axesB;
	const m3d_vec3* incAxes = refIsA ? axesB : axesA;
	const float* refExt = refIsA ? extA : extB;
	const float* incExt = refIsA ? extB : extA;
	m3d_vec3 dRefToInc = refIsA ? d : m3d_neg(d);

	// Reference face normal (world), pointing toward the incident box
	float refSign = m3d_dot(refAxes[refFace], dRefToInc) >= 0.0f ? 1.0f : -1.0f;
	m3d_vec3 refNormal = m3d_scale(refSign, refAxes[refFace]);

	// Incident face: the face of the incident box most anti-parallel to refNormal
	int incFace = 0;
	float minDot = FLT_MAX;
	float incSign = 1.0f;
	for (int i = 0; i < 3; ++i)
	{
		float dot = m3d_dot(incAxes[i], refNormal);
		if (dot < minDot)
		{
			minDot = dot;
			incFace = i;
			incSign = 1.0f;
		}
		if (-dot < minDot)
		{
			minDot = -dot;
			incFace = i;
			incSign = -1.0f;
		}
	}

	// Build the 4 vertices of the incident face in world space
	int iu = (incFace + 1) % 3;
	int iv = (incFace + 2) % 3;
	m3d_vec3 incCenterLocal = m3d_v3_zero();
	((float*)&incCenterLocal)[incFace] = incSign * incExt[incFace];
	m3d_vec3 uL = m3d_v3_zero(), vL = m3d_v3_zero();
	((float*)&uL)[iu] = incExt[iu];
	((float*)&vL)[iv] = incExt[iv];

	ClipVertex poly[8], tmp[8];
	{
		m3d_vec3 cW = m3d_transform_point(xfInc, incCenterLocal);
		m3d_vec3 uW = m3d_mat3_mulv(&rotInc, uL);
		m3d_vec3 vW = m3d_mat3_mulv(&rotInc, vL);
		poly[0] = { m3d_add(m3d_add(cW, uW), vW), 0 };
		poly[1] = { m3d_add(m3d_sub(cW, uW), vW), 1 };
		poly[2] = { m3d_sub(m3d_sub(cW, uW), vW), 2 };
		poly[3] = { m3d_sub(m3d_add(cW, uW), vW), 3 };
	}
	int count = 4;

	// Clip against the 4 side planes of the reference face
	int ru = (refFace + 1) % 3;
	int rv = (refFace + 2) % 3;
	const int sideAxes[2] = { ru, rv };
	for (int k = 0; k < 2; ++k)
	{
		int axisIdx = sideAxes[k];
		m3d_vec3 sn = refAxes[axisIdx];
		float centerProj = m3d_dot(xfRef.p, sn);
		count = ClipPolygon(poly, count, sn, centerProj + refExt[axisIdx], (uint32_t)(2 * k + 4), tmp);
		if (count == 0)
		{
			return;
		}
		count = ClipPolygon(tmp, count, m3d_neg(sn), -(centerProj - refExt[axisIdx]), (uint32_t)(2 * k + 5), poly);
		if (count == 0)
		{
			return;
		}
	}

	// Keep points near or below the reference face
	float faceOffset = m3d_dot(xfRef.p, refNormal) + refExt[refFace];
	ClipVertex kept[8];
	float seps[8];
	int keptCount = 0;
	for (int i = 0; i < count; ++i)
	{
		float sep = m3d_dot(poly[i].v, refNormal) - faceOffset;
		if (sep <= margin)
		{
			kept[keptCount] = poly[i];
			seps[keptCount] = sep;
			++keptCount;
		}
	}
	if (keptCount == 0)
	{
		return;
	}

	// Reduce to at most 4 points: deepest, farthest, then max area both ways
	int order[8];
	int reduced = 0;
	if (keptCount <= 4)
	{
		for (int i = 0; i < keptCount; ++i)
		{
			order[reduced++] = i;
		}
	}
	else
	{
		int i0 = 0;
		for (int i = 1; i < keptCount; ++i)
		{
			if (seps[i] < seps[i0])
			{
				i0 = i;
			}
		}
		int i1 = -1;
		float bestD = -1.0f;
		for (int i = 0; i < keptCount; ++i)
		{
			if (i == i0)
			{
				continue;
			}
			float dd = m3d_length_sq(m3d_sub(kept[i].v, kept[i0].v));
			if (dd > bestD)
			{
				bestD = dd;
				i1 = i;
			}
		}
		int i2 = -1, i3 = -1;
		float bestArea = 0.0f, worstArea = 0.0f;
		for (int i = 0; i < keptCount; ++i)
		{
			if (i == i0 || i == i1)
			{
				continue;
			}
			m3d_vec3 e1 = m3d_sub(kept[i1].v, kept[i0].v);
			m3d_vec3 e2 = m3d_sub(kept[i].v, kept[i0].v);
			float area = m3d_dot(m3d_cross(e1, e2), refNormal);
			if (area > bestArea)
			{
				bestArea = area;
				i2 = i;
			}
			if (area < worstArea)
			{
				worstArea = area;
				i3 = i;
			}
		}
		order[reduced++] = i0;
		order[reduced++] = i1;
		if (i2 >= 0)
		{
			order[reduced++] = i2;
		}
		if (i3 >= 0)
		{
			order[reduced++] = i3;
		}
	}

	// Manifold normal must point A -> B
	m->normal = refIsA ? refNormal : m3d_neg(refNormal);
	uint32_t flipBit = refIsA ? 0u : 0x40000000u;
	for (int k = 0; k < reduced; ++k)
	{
		int i = order[k];
		m3d_vec3 pInc = kept[i].v;
		m3d_vec3 pRef = m3d_mul_add(pInc, -seps[i], refNormal);
		uint32_t id = (kept[i].id & 0x00FFFFFFu) | ((uint32_t)refFace << 24) | flipBit;
		if (refIsA)
		{
			AddManifoldPoint(m, xfA, xfB, pRef, pInc, seps[i], id);
		}
		else
		{
			AddManifoldPoint(m, xfA, xfB, pInc, pRef, seps[i], id);
		}
	}
}

// ---------------------------------------------------------------------------
// Dispatch
// ---------------------------------------------------------------------------

typedef void (*CollideFn)(const Shape&, m3d_transform, const Shape&, m3d_transform, Manifold*);

void CollideShapes(const Shape& shapeA, m3d_transform xfA, const Shape& shapeB, m3d_transform xfB, float margin, Manifold* manifold)
{
	manifold->pointCount = 0;
	manifold->normal = m3d_v3(0.0f, 1.0f, 0.0f);

	if (shapeA.type > shapeB.type)
	{
		// Collide in canonical order, then flip
		Manifold flipped;
		CollideShapes(shapeB, xfB, shapeA, xfA, margin, &flipped);
		manifold->pointCount = flipped.pointCount;
		manifold->normal = m3d_neg(flipped.normal);
		for (int i = 0; i < flipped.pointCount; ++i)
		{
			ManifoldPoint& dst = manifold->points[i];
			const ManifoldPoint& src = flipped.points[i];
			dst = src;
			dst.localAnchorA = src.localAnchorB;
			dst.localAnchorB = src.localAnchorA;
			dst.frictionAnchorA = src.frictionAnchorB;
			dst.frictionAnchorB = src.frictionAnchorA;
		}
		return;
	}

	switch (shapeA.type)
	{
		case M3D_SHAPE_SPHERE:
			switch (shapeB.type)
			{
				case M3D_SHAPE_SPHERE: CollideSpheres(shapeA, xfA, shapeB, xfB, margin, manifold); return;
				case M3D_SHAPE_CAPSULE: CollideSphereCapsule(shapeA, xfA, shapeB, xfB, margin, manifold); return;
				case M3D_SHAPE_BOX: CollideSphereBox(shapeA, xfA, shapeB, xfB, margin, manifold); return;
			}
			return;
		case M3D_SHAPE_CAPSULE:
			switch (shapeB.type)
			{
				case M3D_SHAPE_CAPSULE: CollideCapsules(shapeA, xfA, shapeB, xfB, margin, manifold); return;
				case M3D_SHAPE_BOX: CollideCapsuleBox(shapeA, xfA, shapeB, xfB, margin, manifold); return;
				default: return;
			}
		case M3D_SHAPE_BOX:
			CollideBoxes(shapeA, xfA, shapeB, xfB, margin, manifold);
			return;
	}
}

// ---------------------------------------------------------------------------
// Friction anchor inheritance
// ---------------------------------------------------------------------------

void InheritFrictionAnchors(const Manifold& oldManifold, Manifold& manifold, const Body& a, const Body& b)
{
	constexpr float kSlipReleaseSq = 0.02f * 0.02f;   // stretch beyond this = slip, re-anchor
	constexpr float kFeatureDriftSq = 0.05f * 0.05f;  // sanity: feature must not have jumped

	m3d_vec3 n = manifold.normal;
	for (int32_t i = 0; i < manifold.pointCount; ++i)
	{
		ManifoldPoint& mp = manifold.points[i];
		for (int32_t j = 0; j < oldManifold.pointCount; ++j)
		{
			const ManifoldPoint& old = oldManifold.points[j];
			if (old.id != mp.id)
			{
				continue;
			}
			m3d_vec3 pAold = m3d_add(a.position, m3d_rotate(a.rotation, old.frictionAnchorA));
			m3d_vec3 pBold = m3d_add(b.position, m3d_rotate(b.rotation, old.frictionAnchorB));
			m3d_vec3 pAnew = m3d_add(a.position, m3d_rotate(a.rotation, mp.localAnchorA));

			// stretch of the stored stiction pair, tangentially
			m3d_vec3 gap = m3d_sub(pBold, pAold);
			m3d_vec3 dpt = m3d_sub(gap, m3d_scale(m3d_dot(gap, n), n));
			if (m3d_length_sq(dpt) < kSlipReleaseSq && m3d_length_sq(m3d_sub(pAold, pAnew)) < kFeatureDriftSq)
			{
				mp.frictionAnchorA = old.frictionAnchorA;
				mp.frictionAnchorB = old.frictionAnchorB;
			}
			break;
		}
	}
}

// ---------------------------------------------------------------------------
// Shape ray casts (return fraction in [0, 1] and outward world normal)
// ---------------------------------------------------------------------------

static bool RayCastSphere(m3d_vec3 center, float radius, m3d_vec3 o, m3d_vec3 d, float* fraction, m3d_vec3* normal)
{
	m3d_vec3 s = m3d_sub(o, center);
	float b = m3d_length_sq(s) - radius * radius;
	float rr = m3d_dot(s, d);
	float c = m3d_dot(d, d);
	float sigma = rr * rr - c * b;
	if (sigma < 0.0f || c < 1e-12f)
	{
		return false;
	}
	float t = -(rr + sqrtf(sigma)) / c;
	if (t < 0.0f || t > 1.0f)
	{
		return false;
	}
	*fraction = t;
	*normal = m3d_normalize(m3d_add(s, m3d_scale(t, d)));
	return true;
}

static bool RayCastBoxLocal(m3d_vec3 e, m3d_vec3 o, m3d_vec3 d, float* fraction, m3d_vec3* normalLocal)
{
	float tmin = 0.0f, tmax = 1.0f;
	m3d_vec3 n = m3d_v3_zero();
	const float* op = (const float*)&o;
	const float* dp = (const float*)&d;
	const float* ep = (const float*)&e;
	for (int i = 0; i < 3; ++i)
	{
		if (fabsf(dp[i]) < 1e-12f)
		{
			if (op[i] < -ep[i] || op[i] > ep[i])
			{
				return false;
			}
		}
		else
		{
			float inv = 1.0f / dp[i];
			float t1 = (-ep[i] - op[i]) * inv;
			float t2 = (ep[i] - op[i]) * inv;
			float sign = -1.0f;
			if (t1 > t2)
			{
				float tmpT = t1;
				t1 = t2;
				t2 = tmpT;
				sign = 1.0f;
			}
			if (t1 > tmin)
			{
				tmin = t1;
				n = m3d_v3_zero();
				((float*)&n)[i] = sign;
			}
			tmax = tmax < t2 ? tmax : t2;
			if (tmin > tmax)
			{
				return false;
			}
		}
	}
	if (m3d_length_sq(n) < 0.5f)
	{
		return false; // started inside
	}
	*fraction = tmin;
	*normalLocal = n;
	return true;
}

bool RayCastShape(const Shape& shape, m3d_transform xf, m3d_vec3 origin, m3d_vec3 translation, float* fraction, m3d_vec3* normal)
{
	switch (shape.type)
	{
		case M3D_SHAPE_SPHERE:
			return RayCastSphere(xf.p, shape.radius, origin, translation, fraction, normal);

		case M3D_SHAPE_BOX:
		{
			m3d_vec3 oL = m3d_inv_transform_point(xf, origin);
			m3d_vec3 dL = m3d_inv_rotate(xf.q, translation);
			m3d_vec3 nL;
			if (!RayCastBoxLocal(shape.halfExtents, oL, dL, fraction, &nL))
			{
				return false;
			}
			*normal = m3d_rotate(xf.q, nL);
			return true;
		}

		case M3D_SHAPE_CAPSULE:
		{
			m3d_vec3 oL = m3d_inv_transform_point(xf, origin);
			m3d_vec3 dL = m3d_inv_rotate(xf.q, translation);
			float h = shape.halfHeight;
			float r = shape.radius;
			float best = FLT_MAX;
			m3d_vec3 bestN = m3d_v3_zero();

			// Infinite cylinder x^2 + z^2 = r^2, then clamp |y| <= h
			float aa = dL.x * dL.x + dL.z * dL.z;
			if (aa > 1e-12f)
			{
				float bb = oL.x * dL.x + oL.z * dL.z;
				float cc = oL.x * oL.x + oL.z * oL.z - r * r;
				float sigma = bb * bb - aa * cc;
				if (sigma >= 0.0f)
				{
					float t = (-bb - sqrtf(sigma)) / aa;
					if (t >= 0.0f && t <= 1.0f)
					{
						float y = oL.y + t * dL.y;
						if (fabsf(y) <= h)
						{
							best = t;
							m3d_vec3 hitP = m3d_mul_add(oL, t, dL);
							bestN = m3d_normalize(m3d_v3(hitP.x, 0.0f, hitP.z));
						}
					}
				}
			}

			// Sphere caps
			for (int i = 0; i < 2; ++i)
			{
				m3d_vec3 capCenter = m3d_v3(0.0f, i == 0 ? h : -h, 0.0f);
				float t;
				m3d_vec3 n;
				if (RayCastSphere(capCenter, r, oL, dL, &t, &n) && t < best)
				{
					best = t;
					bestN = n;
				}
			}

			if (best == FLT_MAX)
			{
				return false;
			}
			*fraction = best;
			*normal = m3d_rotate(xf.q, bestN);
			return true;
		}
	}
	return false;
}

} // namespace m3d
