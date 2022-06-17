using Documenter, PGENFiles

makedocs(
    format = Documenter.HTML(),
    sitename = "PGENFiles.jl",
    authors = "Seyoon Ko",
    clean = true,
    debug = true,
    pages = [
        "PGENFiles.jl Tutorial" => "index.md",
        "PGEN format description" => "PGEN_description.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/PGENFiles.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing,
    devbranch = "main"
)
