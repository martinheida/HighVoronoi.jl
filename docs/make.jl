using HighVoronoi
using Documenter

DocMeta.setdocmeta!(HighVoronoi, :DocTestSetup, :(using HighVoronoi); recursive=true)

makedocs(;
    modules=[HighVoronoi],
    authors="Martin Heida",
    repo="https://github.com/martinheida/HighVoronoi.jl/blob/{commit}{path}#{line}",
    sitename="HighVoronoi.jl",
    #versions = ["dev",
    #    "stable",],
    #version = "stable",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://martinheida.github.io/HighVoronoi.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual Voronoi" => Any[
            "Examples Voronoi Generation" => "man/short.md",
            "man/workflowmesh.md",
            "man/geometry.md",
            "man/improving.md",
            "man/boundaries.md",
            "Advanced Options" => "man/advanced.md",
            "man/periodic.md",
            "man/refine.md",
            "Projection operators" => "man/projection.md",
            "man/metapost.md",
            "Sources of errors and loss in performance" => "man/errors.md",
        ],
        "Manual Finite Volume" => Any[
            "Finite Volume Examples" => "man/finitevolumeexample.md",
            "man/workflowfv.md",
            "Finite Volume Tutorial" => "man/finitevolume.md",
            "Functions" => "man/functions.md",
            "man/integrals.md",
            "man/toyfvfile.md",
        ],
        "Intentions of use (EXAMPLES)" => "showcase.md",
    ],
)

deploydocs(;
    repo="github.com/martinheida/HighVoronoi.jl",
    devbranch="main",
    devurl = "v1.1.0",#"v1.0.2",
#    versions = ["stable" => "v^", "v#.#.#",],
)
