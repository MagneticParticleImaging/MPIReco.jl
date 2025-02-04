function MPIReco.reconstruct(name::AbstractString, data::MPIFiles.DMPIFile, cache::Bool = true, modules = [AbstractImageReconstruction, MPIFiles, MPIReco, RegularizedLeastSquares]; kwargs...)
  # Load plan locally
  plan = MPIReco.loadRecoPlan(name, cache, modules; kwargs...)
  # Move to remote
  # TODO Caching
  params = RecoPlan(DaggerReconstructionParameter; worker = MPIFiles.worker(data), algo = plan)
  distr_algo = RecoPlan(DaggerReconstructionAlgorithm; parameter = params)
  # Configure on remote
  setAll!(distr_algo; kwargs...)

  reconstruct(build(distr_algo), data)
end


function MPIReco.reconstruct(algo::AbstractMPIRecoAlgorithm, data::MPIFiles.DMPIFile)
  if myid() != MPIFiles.worker(data)
    @warn "Distributed MPIFile and algorithm exist on different workers. Data has to be transfered over the network"
  else
    data = collect(data.file)
  end
  reconstruct(algo, data)
end