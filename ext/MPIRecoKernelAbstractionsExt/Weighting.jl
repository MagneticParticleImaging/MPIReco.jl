function (params::ProcessResultCache{<:AbstractWeightingParameters})(algo::Type{<:AbstractMPIRecoAlgorithm}, freqs, S, meas, arrType::Type{<:AbstractGPUArray})
  @warn "Caching of weight processing is disabled for GPU processing"
  return params.param(algo, freqs, S, meas, arrType)
end