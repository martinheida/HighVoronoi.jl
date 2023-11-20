# Graphical Output in 2D and 3D

Graphical output can be produced using the `Plots.jl`. This is achieved as follows

## `draw2D`

Example:

```@julia
    VG = VoronoiGeometry(VoronoiNodes(rand(2,5)),cuboid(2,periodic=[]))

    HighVoronoi.draw2D(VG)
```

Optionally, the output can be stored in any format supported by Plots:

```@julia
    HighVoronoi.draw2D(VG,"nice_plot.png")
```


The package provides the following output functions:

```@docs
draw2D
```

You may need the following:
```@docs
PlotBoard
```

## `draw3D`

![Sample 3D plot](./assets/images/3D-Plot.jpg)

The same works with 3 dimensions (you may also pass a customized `PlotBoard`):

```@julia
    VG = VoronoiGeometry(VoronoiNodes(rand(3,5)),cuboid(3,periodic=[]))

    draw3D(VG,"nice_plot_3D.pdf")
```

## Using `PlotlyJS`

```@julia
    using PlotlyJS

    plotly()
    VG = VoronoiGeometry(VoronoiNodes(rand(3,5)),cuboid(3,periodic=[]))

    draw3D(VG)
```

!!! warning " "
    When you use `plotly()` you will not be able to write the output to a file.


## 2D-Output using MetaPost

Similar to LaTeXX, MetaPost is an elegant way to create eps, pdf etc. from a programming language vector graphic code. If you do not have it installed on your PC, you may use the 
MetaPost generator by Troy Henderson: [www.tlhiv.org/mppreview/](http://www.tlhiv.org/mppreview/). However, this link sometimes did not work in the past.


## The MeatPostBoard

These methods are based on the `MetaPostBoard` structure:
```@docs
MetaPostBoard
```

```@docs
MetaPostBoard()
```

