# Advanced Options for controlling the Algorithm and output

## Output options

- `silence=true`: will suppress output by voronoi algorithm and integration, thereby speed up the routine a little bit.
- `printevents=true`: callable only in `VoronoiGeometry(points,...)` this will generate output at the end of the voronoi algorithm with some information on the results and if there where some unexpected events due to non-regular data strucutres.

## Controlling the Voronoi-Algorithm

The optional argument `search_settings::NamedTuple` will help to influence the Voronoi routine. However, if you don't know what for sure what you do, you should not touch these except for `force_irregular_search`. Otherwise don't blame the package if you get strange results or crash the algorithm. 

The following commands are available:
- `recursive` : When set to 'true' this will cause an iteration when verteces are found `out of line`: As long as in some step 'n' verteces are found 
                that primarily belong to vertex `m<n` restart after each full iteration. This is functionality seems to be irrelevant in the latest versions but it may help in case of some very weird data. 
- `variance-tol=1E-20`: when the variance of (distance of a vertex to its nodes)^2 is larger than that value, the vertex candidate will be corrected
- `break_tol=1E-5` : when the afore mentioned variance is even larger than that (actually did not appear in tests on bounded domains so far) this is sign that something goes 
                terribly wrong. Therefore, the vertex is skipped. cases when this happens is a geometry quasi periodic in at least one dimension and unbounded in the others. Typically happens "far away, i.e. 1E200" from the origin.
- `b_nodes_tol=1E-10`: When a vertex has a  distance smaller than that to the boundary, it is considered a boundary vertex. 
- `nodes_tol=1.0E-5,`: When `allow_irregular=true` or `force_irregular_search=true` then this is the threshold below which a node is considered to be part of 
                an irregular vertex.
- `allow_irregular=true`: if `true` the algorithm assumes that irregular verteces may occur and checks for it whenever 
                there is evidence that this might be the case. 
- `force_irregular_search=false`: When `true` the algorithm checks for at every vertex if the vertex is irregular and searches for the irregular nodes. This, however, increases the computational cost. If you assume there will be no irregular verteces, set this field to false
- `correcting=true`: The algorithm checks if the variance of the (renormalized) distance of a vertex to its nodes is larger than `variance_tol` but 
                still smaller than `break_tol`. In this case, it will correct the vertex to reduce the variance as much as possible. Increases the computational cost by less than 1% but is extremely stabilizing for large grids.
