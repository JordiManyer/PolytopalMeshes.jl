
# ===================================================================
# Mesh quality metrics (Section 4.2 of the R3MG paper)
#
# All metrics return per-element values in (0, 1] (higher = better),
# except overlap_factor which is a single scalar for the whole mesh.
# ===================================================================

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

# Cell measures as a collected Vector{Float64}
function _cell_measures(model::DiscreteModel)
  return collect(Float64, CellData.get_cell_measure(Triangulation(model)))
end

# -------------------------------------------------------------------
# Metrics
# -------------------------------------------------------------------

"""
    uniformity_factor(model::DiscreteModel) -> Vector{Float64}

`UF(K) = diam(K) / h` where `h = max_K diam(K)`.

Measures how uniformly sized the elements are.
Values equal to 1 when all elements have the same diameter.
"""
function uniformity_factor(model::DiscreteModel)
  bboxes = _cell_bboxes(model)
  diams  = map(mbr_diameter, bboxes)
  h      = maximum(diams)
  return diams ./ h
end

"""
    circle_ratio(model::DiscreteModel) -> Vector{Float64}

`CR(K) = min_side(MBR(K)) / max_side(MBR(K))`

Aspect-ratio quality measure: 1 for squares/cubes, smaller for elongated
elements.  Bounded in (0, 1].
"""
function circle_ratio(model::DiscreteModel)
  bboxes = _cell_bboxes(model)
  map(bboxes) do b
    max_side = mbr_max_side(b)
    min_side = mbr_min_side(b)
    iszero(max_side) ? 0.0 : min_side / max_side
  end
end

"""
    box_ratio(model::DiscreteModel) -> Vector{Float64}

`BR(K) = |K| / |MBR(K)|`

Ratio of cell measure to its bounding box measure.
Values equal to 1 when the cell exactly fills its bounding box.
"""
function box_ratio(model::DiscreteModel)
  bboxes   = _cell_bboxes(model)
  measures = _cell_measures(model)
  map(zip(measures, bboxes)) do (meas, b)
    box_meas = mbr_measure(b)
    iszero(box_meas) ? 0.0 : meas / box_meas
  end
end

"""
    overlap_factor(model::DiscreteModel) -> Float64

`OF(Ω) = |Ω| / Σ_i |MBR(K_i)|`

Ratio of the total domain measure to the sum of all bounding box
measures.  Values close to 1 mean bounding boxes tightly cover the
domain with little overlap.
"""
function overlap_factor(model::DiscreteModel)
  bboxes       = _cell_bboxes(model)
  domain_meas  = sum(_cell_measures(model))
  sum_box_meas = sum(mbr_measure(b) for b in bboxes)
  iszero(sum_box_meas) ? 0.0 : domain_meas / sum_box_meas
end
