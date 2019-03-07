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
        #"Tracers" => "tracers.md",
        #"Phantoms" => "phantoms.md",
        #"Sequences" => "sequences.md",
        #"Calibration" => "calibration.md",
        #"Datasets" => "datasets.md",
        #"Reconstruction"=> "reconstructions.md"
    ],
)

deploydocs(repo   = "github.com/MagneticParticleImaging/MPIReco.jl.git",
           target = "build")
