# Advanced Options for controlling the Algorithm and output

## Output options

- `silence=true`: will suppress output by voronoi algorithm and integration, thereby speed up the routine a little bit.
- `printevents=true`: callable only in `VoronoiGeometry(points,...)` this will generate output at the end of the voronoi algorithm with some information on the results and if there where some unexpected events due to non-regular data strucutres.

## Controlling the Voronoi-Algorithm

The optional argument `search_settings::NamedTuple` will help to influence the Voronoi routine. However, if you don't know what for sure what you do, you should not touch these. Otherwise don't blame the package if you get strange results or crash the algorithm. Also you should be aware of the fact that parameters once set at generation of a `VoronoiGeometry` are kept for refining, substituting, .... You can, however, locally modify these parameters with `search_settings` in `refine!` or `substitute!`

The following commands are available:
- `variance-tol=1E-20`: when the variance of (distance of a vertex to its nodes)^2 is larger than that value, the vertex candidate will be corrected
- `break_tol=1E-5` : when the afore mentioned variance is even larger than that (actually did not appear in tests on bounded domains so far) this is sign that something goes 
                terribly wrong. Therefore, the vertex is skipped. cases when this happens is a geometry quasi periodic in at least one dimension and unbounded in the others. Typically happens "far away, i.e. 1E200" from the origin.
- `b_nodes_tol=1E-10`: When a vertex has a  distance smaller than that to the boundary, it is considered a boundary vertex. 
