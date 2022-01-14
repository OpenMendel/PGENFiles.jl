using Documenter, PGEN

makedocs(
    format = Documenter.HTML(),
    sitename = "PGEN.jl",
    authors = "Seyoon Ko",
    clean = true,
    debug = true,
    pages = [
        "PGEN.jl Tutorial" => "index.md",
        "PGEN format description" => "PGEN_description.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/PGEN.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    devbranch = "main"
)
