
# ===================================================================
# Generic mesh quality metrics and reporting
#
# All metrics are model-level diagnostics and independent of the
# agglomeration method that produced the mesh.
# ===================================================================

"""
    MeshQualitySummary

Container for quality metrics and aggregate statistics.
"""
struct MeshQualitySummary
  n_cells  :: Int
  uniformity::Vector{Float64}
  circle    ::Vector{Float64}
  box       ::Vector{Float64}
  overlap   ::Float64
end

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

# Cell measures as a collected Vector{Float64}
function _cell_measures(model::DiscreteModel)
  return collect(Float64, CellData.get_cell_measure(Triangulation(model)))
end

@inline function _metric_stats(v::AbstractVector{<:Real})
  return (
    min = minimum(v),
    max = maximum(v),
    mean = Statistics.mean(v),
    median = Statistics.median(v),
    std = Statistics.std(v),
  )
end

function _fmt(x::Real; digits::Int=6)
  return string(round(Float64(x); digits=digits))
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
elements. Bounded in (0, 1].
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
measures. Values close to 1 mean bounding boxes tightly cover the
domain with little overlap.
"""
function overlap_factor(model::DiscreteModel)
  bboxes       = _cell_bboxes(model)
  domain_meas  = sum(_cell_measures(model))
  sum_box_meas = sum(mbr_measure(b) for b in bboxes)
  iszero(sum_box_meas) ? 0.0 : domain_meas / sum_box_meas
end

# -------------------------------------------------------------------
# Driver / report
# -------------------------------------------------------------------

"""
    mesh_quality_summary(model::DiscreteModel) -> MeshQualitySummary

Compute all mesh-quality benchmarks for `model`.
"""
function mesh_quality_summary(model::DiscreteModel)
  uf = uniformity_factor(model)
  cr = circle_ratio(model)
  br = box_ratio(model)
  of = overlap_factor(model)
  return MeshQualitySummary(num_cells(model), uf, cr, br, of)
end

"""
    print_mesh_quality(summary::MeshQualitySummary; io=stdout, title="Mesh Quality")

Print a concise benchmark report.
"""
function print_mesh_quality(summary::MeshQualitySummary;
                            io::IO=stdout,
                            title::AbstractString="Mesh Quality",
                            digits::Int=6)
  uf = _metric_stats(summary.uniformity)
  cr = _metric_stats(summary.circle)
  br = _metric_stats(summary.box)

  println(io, title)
  println(io, "cells: ", summary.n_cells)
  println(io, "uniformity_factor: min=", _fmt(uf.min; digits),
          " max=", _fmt(uf.max; digits),
          " mean=", _fmt(uf.mean; digits),
          " median=", _fmt(uf.median; digits),
          " std=", _fmt(uf.std; digits))
  println(io, "circle_ratio:     min=", _fmt(cr.min; digits),
          " max=", _fmt(cr.max; digits),
          " mean=", _fmt(cr.mean; digits),
          " median=", _fmt(cr.median; digits),
          " std=", _fmt(cr.std; digits))
  println(io, "box_ratio:        min=", _fmt(br.min; digits),
          " max=", _fmt(br.max; digits),
          " mean=", _fmt(br.mean; digits),
          " median=", _fmt(br.median; digits),
          " std=", _fmt(br.std; digits))
  println(io, "overlap_factor:   ", _fmt(summary.overlap; digits))

  return nothing
end

"""
    benchmark_mesh_quality(model::DiscreteModel; io=stdout, title="Mesh Quality")
      -> MeshQualitySummary

Driver utility: compute all quality metrics for `model`, print a report,
and return the computed summary.
"""
function benchmark_mesh_quality(model::DiscreteModel;
                                io::IO=stdout,
                                title::AbstractString="Mesh Quality",
                                digits::Int=6)
  summary = mesh_quality_summary(model)
  print_mesh_quality(summary; io, title, digits)
  return summary
end
