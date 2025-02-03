function MPIReco.reconstruct(name::AbstractString, data::MPIFiles.DMPIFile, cache::Bool = true, modules = [AbstractImageReconstruction, MPIFiles, MPIReco, RegularizedLeastSquares]; kwargs...)
  # Load plan locally
  plan = loadRecoPlan(name, cache, modules; kwargs...)
  # Move to remote
  params = RecoPlan(DaggerReconstructionParameter; worker = data.worker, algo = plan)
  distr_algo = RecoPlan(DaggerReconstructionAlgorithm; parameter = params)
  # Configure on remote
  setAll!(distr_algo; kwargs...)

  reconstruct(build(plan), data)
end


function MPIReco.reconstruct(algo::AbstractMPIRecoAlgorithm, data::MPIFiles.DMPIFile)
  if myid() != data.worker
    @warn "Distributed MPIFile and algorithm exist on different workers"
  end
  reconstruct(algo, collect(data.file))
end