using Documenter, MPIReco

makedocs(
    format = :html,
    modules = [MPIReco],
    sitename = "MPI Reconstruction",
    authors = "Tobias Knopp,...",
    pages = [
        "Home" => "index.md"
        "Usage" => "usage.md"
        #"Scanners" => "scanners.md",
        #"Tracers" => "tracers.md",
        #"Phantoms" => "phantoms.md",
        #"Sequences" => "sequences.md",
        #"Calibration" => "calibration.md",
        #"Datasets" => "datasets.md",
        #"Reconstruction"=> "reconstructions.md"
    ],
)

deploydocs(repo   = "github.com/MagneticParticleImaging/MPIReco.jl.git",
           julia  = "release",
           target = "build",
           deps   = nothing,
           make   = nothing)
