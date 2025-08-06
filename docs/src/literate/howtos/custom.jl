# # Custom Data Processing and Algorithms
include("../../download.jl") #hide
using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf")) #hide
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf")); #hide
# This how-to provides a quick guide on how to customize reconstruction algorithms. For more information on how to extend algorithms and processing steps, we refer to the Explanations section of the documentation
# and to the example section of AbstractImageReconstructions documentation.

# ## Custom Processing Steps
# To implement a custom processing step, we need to add a new parameter struct, in our case we want to extend `AbstractWeightingParameters`.
# As a toy-example, we will implement a weighting strategy in which frequencies are weighting with an alternating sequence of weights:
Base.@kwdef struct AlternatingWeightingParameters <: AbstractWeightingParameters
  alternatingWeights::Vector{Float64}
end
# The `Base.@kwdef` macro generates a keyword-argment constructor with optional default values. Next, we can implement the actual processing function.

# Different algorithms can have different implementations of a given weighting strategy and we can specialise a function on the type of algorithm or an algorithm instance itself.
# The former is helpful for pure functions, i.e. processing steps which solely depend on the given parameter and processing-arguments, not the state of the algorithm.
# The generic weight-process function takes as input the frequencies and the system matrix operator. In our case we only require the frequencies:
function AbstractImageReconstruction.process(::Type{<:AbstractMPIRecoAlgorithm}, params::AlternatingWeightingParameters, freqs::Vector{CartesianIndex{2}}, args...)
  alternatingWeights = Iterators.cycle(params.alternatingWeights) # Infinitely cycle through our weights
  weights = collect(Iterators.take(alternatingWeights, length(freqs)))
  return weights
end

# After implementing our processing function, we can directly use it within our exisintg algorithm blueprints:
params = Dict{Symbol, Any}()
params[:SNRThresh] = 5
params[:frames] = 1:1
params[:minFreq] = 80e3
params[:recChannels] = 1:2
params[:spectralLeakageCorrection] = true
params[:sf] = bSF
params[:reg] = [L2Regularization(0.1f0)]

ourWeighting = AlternatingWeightingParameters(collect(range(0.0, 1.0, length=5)))

# We can use either the high-level interface:
c1 = reconstruct("SinglePatch", b; params..., weightingParams = ourWeighting);

# Or the RecoPlan interface:
plan = MPIRecoPlan("SinglePatch")
setAll!(plan, params)
setAll!(plan, :weightingParams, ourWeighting)
algo = build(plan)
c2 = reconstruct(algo, b)
isapprox(c1.data.data, c2.data.data)

# We can visualize the results of our weighting:
using CairoMakie
fig = heatmap(c1[1, :, :, 1, 1].data.data)
hidedecorations!(fig.axis)
fig

# The `SinglePatch` algorithm which we adapted with our new weighting strategy, actually allows weighting results to be cached.
plan = MPIRecoPlan("SinglePatch")
typeof(plan.parameter.reco.weightingParams)

# By overwritting the weightingParams, we removed the caching layer. If we want to use our adapted reconstruction algoritm with the best performace, we have to retain the caching layer:
plan.parameter.reco.weightingParams.param = ourWeighting
setAll!(plan, params)
algo = build(plan)
c3 = reconstruct(algo, b)
isapprox(c1.data.data, c3.data.data)

# We can then turn out algorithm into a plan and store it somewhere:
plan = toPlan(algo)
AbstractImageReconstruction.clear!(plan)
toTOML(stdout, plan)

# The How-to on [Implement Reconstruction Packages](@ref) shows to make such a plan available to the usual MPIReco interfaces.

# ## Custom Reconstruction Algorithms
# So far, we created new "variant" of an image reconstruction algorithm by adding a new strategy of an existing processing step.
# However, if our new reconstruction algorithm requires overall different data processing, such as for example an X-space reconstruction, we need to implement a new `AbstractMPIRecoAlgorithm`.

# The AbstractImageReconstruction.jl documentation has examples of full algorithm implementations. We recommend reading those and then reading the implementation of the single-patch reconstruction algorithm as a starting point.