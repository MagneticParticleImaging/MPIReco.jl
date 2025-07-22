using MPIReco, Documenter, Literate, DaggerImageReconstruction

# Download data
include("src/download.jl")

# Generate examples
OUTPUT_BASE = joinpath(@__DIR__(), "src", "generated")
INPUT_BASE = joinpath(@__DIR__(), "src", "literate")
for (root, dirs, files) in walkdir(INPUT_BASE)
    for dir in dirs
        OUTPUT = joinpath(OUTPUT_BASE, dir)
        INPUT = joinpath(INPUT_BASE, dir)
        for file in filter(f -> endswith(f, ".jl"), readdir(INPUT))
            Literate.markdown(joinpath(INPUT, file), OUTPUT)
        end
    end
end

makedocs(
    format = Documenter.HTML(prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/MagneticParticleImaging/MPIReco.jl",
        assets=String[],
        collapselevel=1,
    ),
    modules = [MPIReco],
    sitename = "MPIReco",
    authors = "Tobias Knopp,...",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "generated/tutorials/overview.md",
        "Tutorials" => Any[
            "Basic Reconstructions" => "generated/tutorials/basicReconstruction.md",
            "Weighting" => "generated/tutorials/weighting.md",
            #"Background Correction" => "generated/tutorials/background.md",
            "Multi-Contrast" => "generated/tutorials/multiContrast.md",
            "Compression" => "generated/tutorials/matrixCompression.md",
            "Multi-Patch" => "generated/tutorials/multiPatch.md",
            "GPU-Acceleration" => "generated/tutorials/gpuAcceleration.md",
            #"Low-Level" => "generated/tutorials/lowlevel.md",
            "Distributed Reconstruction" => "generated/tutorials/distributed.md"
        ],
        "How to" => Any[
            "Change and Configure Solvers" => "generated/howtos/solvers.md",
            "Enable Caching" => "generated/howtos/caching.md",
            "Implement Custom Data Processing" => "generated/howtos/custom.md",
            "Implement Reconstruction Packages" => "generated/howtos/extensions.md",
        ],
        "Explanations" => Any[
            "Data Structures" => "datastructures.md",
            "MPIRecoPlan" => "generated/explanations/recoplans.md", # What are recoplans, how/where do we load them, how does caching work
            "Imaging Operators" => "generated/explanations/operators.md", # What is required for an operator, how does GPU acceleration work
        ],
        "Reference" => Any[
            "Utility" => "references/utility.md",
            "Parameters" => "references/parameters.md",
            "Single-Patch" => "references/singlepatch.md",
            "Multi-Patch" => "references/multipatch.md",
        ]
    ],
    warnonly = [:missing_docs]
)

deploydocs(repo   = "github.com/MagneticParticleImaging/MPIReco.jl.git",
           target = "build")
