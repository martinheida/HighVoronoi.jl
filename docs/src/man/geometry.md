# Creating and handling the Voronoi geometry

The most basic thing is the creation of a list of Points. We advise to use the following:

```@docs
VoronoiNodes(x::Matrix)
```

The creation and storage of Voronoi geometry data is handled by the following class. 

```@docs
VoronoiGeometry{T}
```

To create a Voronoi mesh it is most convenient to call either of the following methods

```@docs
VoronoiGeometry()
```

## Storage

```@docs
write_jld()
```

```@docs
load_Voronoi_info()
```

## Extraction of `VoronoiData` data for further processing

```@docs
VoronoiData
```

