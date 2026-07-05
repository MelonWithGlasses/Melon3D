// Sorted spatial hash grid broadphase.
//
// Every step the grid is rebuilt from scratch: each body's expanded AABB is
// binned into the cells it overlaps and the (cellKey, body) entries are
// sorted. Candidate pairs come from bodies sharing a cell; pairs spanning
// several cells are deduplicated with the min-corner-of-intersection rule.
// Bodies covering too many cells (e.g. the ground) go to a coarse list that
// is tested against everything - O(oversized * n), cheap for a handful of
// large statics.
#include "m3d_internal.h"

#include <algorithm>
#include <math.h>

namespace m3d
{

static inline int64_t CellCoord(float v, float invCell)
{
	return (int64_t)floorf(v * invCell);
}

static inline uint64_t PackCell(int64_t x, int64_t y, int64_t z)
{
	const int64_t offset = 1 << 20; // 21 bits per axis
	uint64_t ux = (uint64_t)(x + offset) & 0x1FFFFF;
	uint64_t uy = (uint64_t)(y + offset) & 0x1FFFFF;
	uint64_t uz = (uint64_t)(z + offset) & 0x1FFFFF;
	return (ux << 42) | (uy << 21) | uz;
}

uint64_t SpatialGrid::CellKeyOf(m3d_vec3 p) const
{
	float inv = 1.0f / m_cellSize;
	return PackCell(CellCoord(p.x, inv), CellCoord(p.y, inv), CellCoord(p.z, inv));
}

void SpatialGrid::BuildTier(float cellSize, const std::vector<m3d_aabb>& aabbs, const std::vector<uint8_t>& include,
							std::vector<Entry>& entries, std::vector<int32_t>& oversized, std::vector<bool>& inGrid)
{
	entries.clear();
	oversized.clear();
	inGrid.assign(aabbs.size(), false);

	float inv = 1.0f / cellSize;
	for (int32_t i = 0; i < (int32_t)aabbs.size(); ++i)
	{
		if (!include[i])
		{
			continue;
		}
		const m3d_aabb& box = aabbs[i];
		int64_t x0 = CellCoord(box.lowerBound.x, inv);
		int64_t y0 = CellCoord(box.lowerBound.y, inv);
		int64_t z0 = CellCoord(box.lowerBound.z, inv);
		int64_t x1 = CellCoord(box.upperBound.x, inv);
		int64_t y1 = CellCoord(box.upperBound.y, inv);
		int64_t z1 = CellCoord(box.upperBound.z, inv);

		int64_t span = (x1 - x0 + 1) * (y1 - y0 + 1) * (z1 - z0 + 1);
		if (span > kOversizedCellSpan)
		{
			oversized.push_back(i);
			continue;
		}
		inGrid[i] = true;
		for (int64_t x = x0; x <= x1; ++x)
		{
			for (int64_t y = y0; y <= y1; ++y)
			{
				for (int64_t z = z0; z <= z1; ++z)
				{
					entries.push_back({ PackCell(x, y, z), i });
				}
			}
		}
	}

	std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
		return a.key < b.key || (a.key == b.key && a.body < b.body);
	});
}

void SpatialGrid::BuildActive(float cellSize, const std::vector<m3d_aabb>& aabbs, const std::vector<uint8_t>& isActive)
{
	m_cellSize = cellSize;
	BuildTier(cellSize, aabbs, isActive, m_active, m_activeOversized, m_inActive);
}

void SpatialGrid::BuildStable(float cellSize, const std::vector<m3d_aabb>& aabbs, const std::vector<uint8_t>& isStable)
{
	m_cellSize = cellSize;
	BuildTier(cellSize, aabbs, isStable, m_stable, m_stableOversized, m_inStable);
}

m3d_aabb ComputeShapeAABB(const Shape& shape, m3d_transform xf)
{
	switch (shape.type)
	{
		case M3D_SHAPE_SPHERE:
		{
			m3d_vec3 r = m3d_v3(shape.radius, shape.radius, shape.radius);
			return { m3d_sub(xf.p, r), m3d_add(xf.p, r) };
		}
		case M3D_SHAPE_CAPSULE:
		{
			m3d_vec3 axis = m3d_rotate(xf.q, m3d_v3(0.0f, shape.halfHeight, 0.0f));
			m3d_vec3 p1 = m3d_add(xf.p, axis);
			m3d_vec3 p2 = m3d_sub(xf.p, axis);
			m3d_vec3 r = m3d_v3(shape.radius, shape.radius, shape.radius);
			return { m3d_sub(m3d_min(p1, p2), r), m3d_add(m3d_max(p1, p2), r) };
		}
		case M3D_SHAPE_BOX:
		default:
		{
			m3d_mat3 rot = m3d_quat_to_mat3(xf.q);
			m3d_vec3 e = shape.halfExtents;
			m3d_vec3 ext;
			ext.x = fabsf(rot.cx.x) * e.x + fabsf(rot.cy.x) * e.y + fabsf(rot.cz.x) * e.z;
			ext.y = fabsf(rot.cx.y) * e.x + fabsf(rot.cy.y) * e.y + fabsf(rot.cz.y) * e.z;
			ext.z = fabsf(rot.cx.z) * e.x + fabsf(rot.cy.z) * e.y + fabsf(rot.cz.z) * e.z;
			return { m3d_sub(xf.p, ext), m3d_add(xf.p, ext) };
		}
	}
}

} // namespace m3d
