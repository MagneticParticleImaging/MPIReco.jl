# # Implement Reconstruction Packages

# In the [Custom Data Processing and Algorithms](@ref) how-to, we showed how to create custom parameters, processing steps and algorithms.
# We also showed that we can create empty `RecoPlans` with new custom structures. The only thing missing to fully make such blueprints available
# to the usual MPIReco interface is to essentially tell MPIReco about our custom structures and blueprints.

# Let's implement a small Julia package for our weighting strategy:
module CustomWeighting
  using MPIReco
  using MPIReco.AbstractImageReconstruction

  export AlternatingWeightingParameters
  Base.@kwdef struct AlternatingWeightingParameters <: AbstractWeightingParameters
    alternatingWeights::Vector{Float64}
  end

  function AbstractImageReconstruction.process(::Type{<:AbstractMPIRecoAlgorithm}, params::AlternatingWeightingParameters, freqs::Vector{CartesianIndex{2}}, args...)
    alternatingWeights = Iterators.cycle(params.alternatingWeights) # Infinitely cycle through our weights
    weights = collect(Iterators.take(alternatingWeights, length(freqs)))
    return weights
  end
end

# To use this module, a user could simply import both MPIReco and CustomWeighting:
using MPIReco, .CustomWeighting
plan = MPIRecoPlan("SinglePatch")
setAll!(plan, :weightingParams, AlternatingWeightingParameters(collect(range(0.0, 1.0, length=5))))

# However, they won't be able to do something like
# ```julia
# reconstruct("SinglePatchAlternating", b; kwargs...)
# ```
# just yet.

# We first have to store a blueprint in some directory `dir`. This could be an directory within our `CustomWeighting` package or an expected folder in the users filesystem. Secondly, we need to tell MPIReco which modules
# it should consider when loading a blueprint. We can do both these things in the initialization function of our package:
# ```julia
# module CustomWeighting
#   # ... same code as above
# 
#   dir = joinpath(@__DIR__(), "..", "Plans")
#   function __init__()
#     addRecoPlanPath(dir)
#     addRecoPlanModule(CustomWeighting)
#   end
# end
# ```

# Afterwards a user can simply invoke our reconstruction algorithms just like they would the provided reconstruction algorithms.