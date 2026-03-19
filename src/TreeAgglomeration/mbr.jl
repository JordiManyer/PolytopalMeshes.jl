
# ===================================================================
# MBR — Minimum Bounding Rectangle (axis-aligned bounding box)
# ===================================================================

"""
    BBox{D,T}

Axis-aligned bounding box in `D` dimensions with coordinate type `T`.
Stored as a pair of corner tuples `(lo, hi)`.
"""
struct BBox{D, T <: AbstractFloat}
  lo :: NTuple{D, T}
  hi :: NTuple{D, T}
end

BBox(lo::Point{D,T}, hi::Point{D,T}) where {D,T} = BBox{D,T}(lo.data, hi.data)

# Compute the MBR from an array of mesh-cell vertex coordinates
function BBox(cell_coords)
  plo, phi = ReferenceFEs.get_bounding_box(cell_coords)
  return BBox(plo, phi)
end

# ---- Geometric operations ----------------------------------------

"""Return the combined MBR that encloses both `a` and `b`."""
@inline function combine(a::BBox{D,T}, b::BBox{D,T}) where {D,T}
  BBox{D,T}(
    ntuple(d -> min(a.lo[d], b.lo[d]), Val(D)),
    ntuple(d -> max(a.hi[d], b.hi[d]), Val(D)),
  )
end

"""D-dimensional measure (length / area / volume) of `b`."""
@inline function mbr_measure(b::BBox{D,T}) where {D,T}
  prod(ntuple(d -> b.hi[d] - b.lo[d], Val(D)))
end

"""Length of the main diagonal of `b`."""
@inline function mbr_diameter(b::BBox{D,T}) where {D,T}
  sqrt(sum(ntuple(d -> (b.hi[d] - b.lo[d])^2, Val(D))))
end

"""Shortest side length of `b`."""
@inline function mbr_min_side(b::BBox{D,T}) where {D,T}
  minimum(ntuple(d -> b.hi[d] - b.lo[d], Val(D)))
end

"""Longest side length of `b`."""
@inline function mbr_max_side(b::BBox{D,T}) where {D,T}
  maximum(ntuple(d -> b.hi[d] - b.lo[d], Val(D)))
end

"""Centroid of `b` as an NTuple."""
@inline function mbr_center(b::BBox{D,T}) where {D,T}
  ntuple(d -> (b.lo[d] + b.hi[d]) / 2, Val(D))
end

"""D-dimensional overlap volume between `a` and `b` (0 if they don't intersect)."""
@inline function mbr_overlap(a::BBox{D,T}, b::BBox{D,T}) where {D,T}
  prod(ntuple(d -> max(zero(T), min(a.hi[d], b.hi[d]) - max(a.lo[d], b.lo[d])), Val(D)))
end
