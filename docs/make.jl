push!(LOAD_PATH,".../src/")
using Revise, Documenter, HighVoronoi

makedocs(sitename="HighVoronoi Documentation",
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
#        prettyurls = !("local" in ARGS),
#        canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
#        assets = ["assets/favicon.ico"],
#        analytics = "UA-136089579-2",
#        highlights = ["yaml"],
        ansicolor = true,
    ),
    authors = "Martin Heida.",
    pages = [
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
#        "Library" => Any[
#            "Public" => "lib/public.md",
#            "Internals" => map(
#                s -> "lib/internals/$(s)",
#                sort(readdir(joinpath(@__DIR__, "src/lib/internals")))
#            ),
#        ],
#        "contributing.md",
    ],
)

#=
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Guide" => "man/guide.md",
            "man/examples.md",
            "man/syntax.md",
            "man/doctests.md",
            "man/latex.md",
            hide("man/hosting.md", [
                "man/hosting/walkthrough.md"
            ]),
            "man/other-formats.md",
        ],
        "showcase.md",
        "Library" => Any[
            "Public" => "lib/public.md",
            "Internals" => map(
                s -> "lib/internals/$(s)",
                sort(readdir(joinpath(@__DIR__, "src/lib/internals")))
            ),
        ],
        "contributing.md",
    ],

    =#