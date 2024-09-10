# Voronoi: Raycast methods

During the develop of `HighVoronoi.jl` I developed several `Raycast`-methods. These methods are central to the Voronoi algorithm.

## [Classical Method](@id classicraycast)

This is the most classical method introduced in [PREPRINT](http://www.wias-berlin.de/preprint/3041/wias_preprints_3041.pdf). It basically works only on generators in general position (when a vertex is created by exactly $d+1$ generators) and can be called using `method=RCOriginal` in the `search_settings`.

## [Statistically Boosted Method](@id boostedraycast)

This method is called with `method=RCNonGeneral` and speeds up the classical method in the following way: It performs two classical steps as in [PREPRINT](http://www.wias-berlin.de/preprint/3041/wias_preprints_3041.pdf) and then performs an `inrange` search and finally sorts out all non-generators. This provides a significant boost when the generators are in non-general position (more than $d+1$ generators for one vertex)

## [Nested Method](@id nestedraycast)

This method is called with `method=RCCombined` and performs a deep hacking of the nearest-neighbor search, i.e. it made it necessary to include a modified version of `NearestNeighbors.jl` into the source code of `HighVoronoi.jl`. The result is a nearest neighbor search that basically converges the raycast method within one single search. It is very fast on bounded and periodic domains but may take much longer on unbounded domains at the periphery of the point set.


## Suggested usage

It is highly recommended to use `method=RCCombined` (the standard setting) and switch to `method=RCNonGeneral` in pathological situation respectively `method=RCNonGeneral` only if the generators are "nicely" distributed.
