# # Enable Caching
include("../../download.jl") #hide

# Image reconstructions implemented with AbstractImageReconstruction.jl are composed of several individual processing steps which form the whole (computationally expensive) image recosntruction procses.
# Often time one wants to slightly modify reconstruction parameters, for example when searching for good regularization parameters.
# AbstractImageReconstruction.jl offers a caching option in which algorithms can reuse intermediate processing results, as long as those are unaffected by the parameter changes.

# This caching mechanism has to be explicitly added in the implementation of an algorithm and an algorithms RecoPlan. All algorithms constructed from the same plan, can access the same cache.
# To make this behaviour available to the high-level reconstruct interface, MPIReco can cache RecoPlans between reconstructions:
using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))
params = Dict{Symbol, Any}()
params[:SNRThresh] = 5
params[:frames] = 1:1
params[:minFreq] = 80e3
params[:recChannels] = 1:2
params[:spectralLeakageCorrection] = true
params[:sf] = bSF
params[:reg] = [L2Regularization(0.1f0)];

# Caching is disabled by default and if you modify the blueprint of a plan, for example by providing a different weighting strategy.
# To opt into caching you call a reconstruction with a true flag as the third argument:
initial = @elapsed reconstruct("SinglePatch", b, true; params...)
second = @elapsed reconstruct("SinglePatch", b, true; params...)
initial/second # speedup

# Changing a parameter which affects the cache result, in this case the loading of the system matrix, invalidates the cache:
third = @elapsed reconstruct("SinglePatch", b, true; params..., SNRThresh = 2)
initial/third

# It is possible to manually empty the MPIRecos cache of blueprints and thus release any large cached results:
emptyRecoCache!();
# The number of cached reconstruction plans can be set via the `MPIRECO_CACHE_SIZE` environment variable. Note that restarting is necessary for it to take any effect.

# When directly constructing a plan, all algorithms build from it can benefit from the plans caching:
plan = MPIRecoPlan("SinglePatch")
setAll!(plan, params)
results = []
for λ in [0.1f0, 0.5f0, 1.0f0]
  setAll!(plan, :reg, [L2Regularization(λ)])
  c = reconstruct(build(plan), b)
  push!(results, c[1, :, :, 1, 1].data)
end
using CairoMakie #hide
fig = Figure();
hidedecorations!(heatmap(fig[1, 1], results[1], axis = (title = "λ = 0.1",)).axis)
hidedecorations!(heatmap(fig[1, 2], results[2], axis = (title = "λ = 0.5",)).axis)
hidedecorations!(heatmap(fig[1, 3], results[3], axis = (title = "λ = 1.0",)).axis)
fig
# This type of caching bypasses MPIRecos cache and results are only freed if there is no reference to the plan available anymore:
plan = nothing