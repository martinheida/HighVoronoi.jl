using HighVoronoi
using Documenter

DocMeta.setdocmeta!(HighVoronoi, :DocTestSetup, :(using HighVoronoi); recursive=true)

makedocs(;
    modules=[HighVoronoi],
    authors="Martin Heida",
    repo="https://github.com/martinheida/HighVoronoi.jl/blob/{commit}{path}#{line}",
    sitename="HighVoronoi.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://martinheida.github.io/HighVoronoi.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "In short..." => "man/short.md",
            "man/geometry.md",
            "man/boundaries.md",
            "man/refine.md",
            "man/periodic.md",
            "FunctionComposer" => "man/functions.md",
            "Finite Volume Problems" => "man/finitevolume.md",
            "Projection operators" => "man/projection.md",
            "man/metapost.md",
            "Sources of errors and loss in performance" => "man/errors.md",
        ],
        "Intentions of use (EXAMPLES)" => "showcase.md",
    ],
)

deploydocs(;
    repo="github.com/martinheida/HighVoronoi.jl",
    devbranch="main",
)
