using Documenter, MPIReco

makedocs(
    format = :html,
    modules = [MPIReco],
    sitename = "MPI Reconstruction",
    authors = "Tobias Knopp,...",
    pages = [
        "Home" => "index.md",
        "Installation" => "usage.md",
        "Overview" => "overview.md",
        "Multi-patch" => "multiPatch.md"
    ],
    html_prettyurls = false, #!("local" in ARGS),
)

deploydocs(repo   = "github.com/MagneticParticleImaging/MPIReco.jl.git",
           target = "build")
