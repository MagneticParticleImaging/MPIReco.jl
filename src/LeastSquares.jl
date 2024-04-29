export LeastSquaresParameters
# TODO this could be moved to AbstractImageReconstruction, depends on how MRIReco.jl structures its data arrays
abstract type AbstractSolverParameters <: AbstractMPIRecoParameters end

export LeastSquaresParameters
Base.@kwdef struct LeastSquaresParameters{L<:AbstractLinearSolver, O, M, R<:AbstractRegularization, P<:AbstractSolverParameters} <: AbstractMPIRecoParameters
  solver::Type{L} = Kaczmarz
  op::O = nothing
  S::M
  reg::Vector{R} 
  solverParams::P
end

# TODO place weights and more

export SimpleSolverParameters
Base.@kwdef struct SimpleSolverParameters <: AbstractSolverParameters
  iterations::Int64=10
  enforceReal::Bool=true
  enforcePositive::Bool=true
  normalizeReg::AbstractRegularizationNormalization = SystemMatrixBasedNormalization()
end
export ConstraintMaskedSolverParameters
Base.@kwdef struct ConstraintMaskedSolverParameters{P<:AbstractSolverParameters} <: AbstractSolverParameters
  constraintMask::Vector{Bool}
  params::P
end

function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::LeastSquaresParameters, u::Array)

  N = size(params.S, 2)
  M = div(length(params.S), N)
  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  reg, args = prepareRegularization(params.reg, params)
  args[:reg] = reg

  solv = createLinearSolver(params.solver, params.S; args...)

  for l=1:L
    d = solve!(solv, u[:, l])
    if !isnothing(params.op)
      d[:] = params.op*d
    end
    c[:, l] = real(d)
  end

  return c

end

function prepareRegularization(reg::Vector{R}, regLS::LeastSquaresParameters) where R<:AbstractRegularization
  args = toKwargs(regLS.solverParams)
  params = regLS.solverParams

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