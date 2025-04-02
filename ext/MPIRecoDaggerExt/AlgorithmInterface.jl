"""
    reconstruct(name::AbstractString, data::MPIFiles.DMPIFile, cache::Bool = false, modules = [AbstractImageReconstruction, MPIFiles, MPIReco, RegularizedLeastSquares]; kwargs...)

Perform a reconstruction on the worker specified by `DMPIFile`. Cache entries are local to the respective workers.
"""
function MPIReco.reconstruct(name::AbstractString, data::MPIFiles.DMPIFile, cache::Bool = false, modules = [AbstractImageReconstruction, MPIFiles, MPIReco, RegularizedLeastSquares]; kwargs...)
  planfile = MPIReco.planpath(name)

  # If the user disables caching or changes the plan structure we bypass the cache
  kwargValues = values(values(kwargs))
  if !cache || any(val -> isa(val, AbstractRecoPlan) || isa(val, AbstractImageReconstructionParameters), kwargValues)
    plan = loadDaggerPlan(planfile, modules; worker = MPIFiles.worker(data))
  else
    key = hash(planfile, hash(mtime(planfile), hash(MPIFiles.worker(data))))
    plan = get!(MPIReco.recoPlans, key) do
      loadDaggerPlan(planfile, modules; worker = MPIFiles.worker(data))
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