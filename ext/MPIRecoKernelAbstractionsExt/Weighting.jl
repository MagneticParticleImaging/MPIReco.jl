function AbstractImageReconstruction.process(algo::Type{<:AbstractMPIRecoAlgorithm}, params::ProcessResultCache{<:AbstractWeightingParameters}, freqs, sf, S, meas, arrType::Type{<:AbstractGPUArray})
  @warn "Caching of weight processing is disabled for GPU processing"
  return process(algo, params.param, freqs, sf, S, meas, arrType)
end