module MPIRecoDaggerExt

using MPIReco, DaggerImageReconstruction, Dagger, DaggerImageReconstruction.Distributed
using MPIFiles

const MPIFilesDaggerExt = Base.get_extension(MPIFiles, :MPIFilesDaggerExt)

if isnothing(MPIFilesDaggerExt)
  error("Could not retrieve MPIFiles extension for Dagger.jl")
end

include("AlgorithmInterface.jl")

end