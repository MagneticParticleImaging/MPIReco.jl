module MPIRecoDaggerExt

using MPIReco, DaggerImageReconstruction, Dagger, DaggerImageReconstruction.Distributed
using MPIReco.AbstractImageReconstruction
using MPIReco.MPIFiles

#https://github.com/JuliaLang/julia/pull/56368 didn't work
#const MPIFilesDaggerExt = Base.get_extension(MPIFiles, :MPIFilesDaggerExt)
#if isnothing(MPIFilesDaggerExt)
#  error("Could not retrieve MPIFiles extension for Dagger.jl")
#end

if isdefined(MPIFiles, :DMPIFile)
  include("AlgorithmInterface.jl")
  include("DaggerRecoPlan.jl")
else
  @warn "MPIFiles version does not export `DMPIFile`s, automatic distributed reconstruction is not possible"
end

end