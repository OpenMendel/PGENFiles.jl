using Documenter, PGEN

makedocs(
    format = Documenter.HTML(),
    sitename = "PGEN.jl",
    authors = "Seyoon Ko",
    clean = true,
    debug = true,
    pages = [
        "index.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/PGEN.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    devbranch = "main"
)