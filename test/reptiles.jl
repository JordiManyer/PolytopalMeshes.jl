
using Gridap
using Gridap.ReferenceFEs, Gridap.Adaptivity, Gridap.Geometry, Gridap.Arrays
using Gridap.Fields, Gridap.TensorValues

using DrWatson

using GridapBDDC
using TikzPictures

# L refinement

rr = LReptileRefinementRule()
TikzPicture(rr.ref_grid.grid_topology;options="color=white, scale=5",label_offset=0.02)

cmodel = rr.ref_grid
fmodel = reptile_refinement(cmodel, rr)
for i in 1:3
  fmodel = reptile_refinement(fmodel.model, rr)
end
writevtk(fmodel.model, datadir("l_reptile_model"))

# Sphinx refinement

rr = SphinxReptileRefinementRule()
TikzPicture(rr.ref_grid.grid_topology;options="color=white, scale=5",label_offset=0.02)
