function MPIReco.reconstruct(name::AbstractString, data::MPIFiles.DMPIFile, cache::Bool = true, modules = [AbstractImageReconstruction, MPIFiles, MPIReco, RegularizedLeastSquares]; kwargs...)
  planfile = AbstractImageReconstruction.planpath(MPIReco, name)

  # If the user disables caching or changes the plan structure we bypass the cache
  kwargValues = values(values(kwargs))
  if !cache || any(val -> isa(val, AbstractRecoPlan) || isa(val, AbstractImageReconstructionParameters), kwargValues)
    plan = distribute_plan(MPIReco.loadRecoPlan(planfile, modules), MPIFiles.worker(data))
  else
    key = hash(planfile, hash(mtime(planfile), hash(MPIFiles.worker(data))))
    plan = get!(MPIReco.recoPlans, key) do
      distribute_plan(MPIReco.loadRecoPlan(planfile, modules), MPIFiles.worker(data))
    end
  end
  # Configure on remote
  setAll!(plan; kwargs...)

  reconstruct(build(plan), data)
end

function distribute_plan(plan, worker)
  params = RecoPlan(DaggerReconstructionParameter; worker = worker, algo = plan)
  distr_algo = RecoPlan(DaggerReconstructionAlgorithm; parameter = params)
  return distr_algo
end


function MPIReco.reconstruct(algo::AbstractMPIRecoAlgorithm, data::MPIFiles.DMPIFile)
  if myid() != MPIFiles.worker(data)
    @warn "Distributed MPIFile and algorithm exist on different workers. Data has to be transfered over the network"
  else
    data = collect(data.file)
  end
  reconstruct(algo, data)
end