
# [Convex Hull](@id convexhull)

Use the following to compute the convex hull of the `VoronoiNodes` `xs`:

```@julia
cv = ConvexHull(xs,"computing convex hull: "; nthreads=Threads.nthreads(), method = RCOriginal)
```

- Note that only the following methods are available: `RCOriginal`, `RCNonGeneral`. Again, `RCOriginal` is faster and for nodes in general position (that means there are only $d$ points in $\mathbb{R}^d$ that generate one surface element) while  `RCNonGeneral` can be used for nodes in non-general position.
- `nthreads` will be initialized with `Threads.nthreads()` by default. If you want less threads, provide a positive `Int64`. You can set the string to ` "" ` in order to have only a progress bar. 
- `length(cv)` will give you the number of surface elements
- `cv[i]` will then give you `(sig,r,u)` where `sig` is a list of generating nodes, `r` is a point on the plane equidistant to all `sig` and `u` is the outer normal of that plane. 
