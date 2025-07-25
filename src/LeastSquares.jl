export LeastSquaresParameters
# TODO this could be moved to AbstractImageReconstruction, depends on how MRIReco.jl structures its data arrays
abstract type AbstractSolverParameters{AbstractLinearSolver} <: AbstractMPIRecoParameters end

export LeastSquaresParameters
Base.@kwdef struct LeastSquaresParameters{L<:AbstractLinearSolver, O, R<:AbstractRegularization, P<:AbstractSolverParameters{L}, W} <: AbstractMPIRecoParameters
  op::O = nothing
  reg::Vector{R} 
  solverParams::P
  weights::W = nothing
end

# TODO place weights and more

export SimpleSolverParameters
Base.@kwdef struct SimpleSolverParameters <: AbstractSolverParameters{Kaczmarz}
  iterations::Int64=10
  enforceReal::Bool=true
  enforcePositive::Bool=true
  normalizeReg::AbstractRegularizationNormalization = SystemMatrixBasedNormalization()
end
export ConstraintMaskedSolverParameters
Base.@kwdef struct ConstraintMaskedSolverParameters{S, P<:AbstractSolverParameters{S}} <: AbstractSolverParameters{S}
  constraintMask::Vector{Bool}
  params::P
end
export ElaborateSolverParameters
Base.@kwdef mutable struct ElaborateSolverParameters{SL} <: AbstractSolverParameters{SL}
  solver::Type{SL} = Kaczmarz
  iterations::Int64 = 10
  enforceReal::Bool = true
  enforcePositive::Bool = true
  kwargWarning::Bool = true
  # Union of all kwargs
  normalizeReg::AbstractRegularizationNormalization = SystemMatrixBasedNormalization()
  randomized::Union{Nothing, Bool, Symbol} = nothing
  greedy_randomized::Union{Nothing, Bool} = nothing
  subMatrixFraction::Union{Nothing, Float64} = nothing
  seed::Union{Nothing, Int64} = nothing
  shuffleRows::Union{Nothing, Bool} = false
  rho::Union{Nothing, Float64} = nothing
  vary_rho::Union{Nothing, Symbol} = nothing
  theta::Union{Nothing, Float64} = nothing
  restart::Union{Nothing, Symbol} = nothing
  regTrafo::Union{Nothing, Vector{Union{AbstractArray, AbstractLinearOperator}}} = nothing
  relTol::Union{Nothing, Float64} = nothing
  absTol::Union{Nothing, Float64} = nothing
  tolInner::Union{Nothing, Float64} = nothing
  iterationsCG::Union{Nothing, Int64} = nothing
  iterationsInner::Union{Nothing, Int64} = nothing
end
Base.propertynames(params::ElaborateSolverParameters{SL}) where SL = union([:solver, :iterations, :enforceReal, :enforcePositive], getSolverKwargs(SL))
Base.propertynames(params::RecoPlan{ElaborateSolverParameters}) = union([:solver, :iterations, :enforceReal, :enforcePositive], ismissing(params.solver) ? getSolverKwargs(Kaczmarz) : getSolverKwargs(params.solver))

getSolverKwargs(::Type{SL}) where SL <: AbstractLinearSolver = intersect(union(Base.kwarg_decl.(methods(SL))...), fieldnames(ElaborateSolverParameters))

function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::LeastSquaresParameters{SL}, S, u::AbstractArray) where SL

  N = size(S, 2)
  M = div(length(S), N)
  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  reg, args = prepareRegularization(params.reg, params)
  args[:reg] = reg

  if !isnothing(params.weights)
    S = ProdOp(WeightingOp(params.weights), S)
    u = params.weights.*u
  end
  SHS = prepareNormalSF(SL, S)
  args[:AHA] = SHS

  solv = createLinearSolver(SL, S; filter(entry -> !isnothing(entry.second), args)...)

  for l=1:L
    d = solve!(solv, u[:, l])
    if !isnothing(params.op)
      d[:] = params.op*d
    end
    c[:, l] = Array(real(d))
  end

  return c

end

function prepareRegularization(reg::Vector{R}, regLS::LeastSquaresParameters) where R<:AbstractRegularization
  params = regLS.solverParams
  args = toKwargs(params)
  if haskey(args, :solver)
    pop!(args, :solver)
  end

  result = AbstractRegularization[]
  push!(result, reg...)

  # Add further constraints
  if params.enforceReal && params.enforcePositive
    if !any([item isa PositiveRegularization for item ∈ reg])
      push!(result, PositiveRegularization())
    end

    filter!(x -> !(x isa RealRegularization), reg)
  elseif params.enforceReal
    if !any([item isa RealRegularization for item ∈ reg])
      push!(result, RealRegularization())
    end
  end

  pop!(args, :enforceReal)
  pop!(args, :enforcePositive)

  # Add sparsity op
  if !isnothing(regLS.op)
    result = map(r -> TransformedRegularization(r, regLS.op), result)
  end

  return result, args
end
#=
function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::LeastSquaresParameters, threadInput::MultiThreadedInput)
  
  scheduler = threadInput.scheduler
  data = threadInput.inputs
  u = threadInput.inputs[1]
  
  N = size(params.S, 2)
  M = div(length(params.S), N)
  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  args = toKwargs(params.solverParams)
  reg = prepareRegularization(params.reg, params)
  args[:reg] = reg
  args[:sparseTrafo] = params.op
  solv = createLinearSolver(params.solver, params.S; args..., weights = params.weights)

  # Distribute frames onto nthreads tasks
  numThreads = AbstractImageReconstruction.nthreads(scheduler)
  threadFrames = Vector{UnitRange{Int64}}()
  solvers = Vector{AbstractLinearSolver}()
  for block in Iterators.partition(1:L, ceil(Int64, L/numThreads))
    push!(threadFrames, block)
    push!(solvers, copy(solv))
  end

  function threadSolve(solver, frames)
    for frame in frames
      d = solve!(solver, u[:, frame])
      if !isnothing(params.op)
        d[:] = params.op*d
      end
      c[:,frame] = real(d)
    end
    return nothing
  end


  for (i, frames) in enumerate(threadFrames)
    put!(scheduler, threadSolve, solvers[i], frames)
  end
  result = nothing
  for frames in threadFrames
    take!(scheduler)
  end
  return c
end
=#