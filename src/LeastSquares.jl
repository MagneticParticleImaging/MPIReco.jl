export LeastSquaresParameters
# TODO this could be moved to AbstractImageReconstruction, depends on how MRIReco.jl structures its data arrays
abstract type AbstractSolverParameters <: AbstractMPIRecoParameters end

export LeastSquaresParameters
Base.@kwdef struct LeastSquaresParameters{L<:AbstractLinearSolver, O, M, R<:AbstractRegularization, P<:AbstractSolverParameters} <: AbstractMPIRecoParameters
  solver::Type{L}
  op::O
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

function process(t::Type{<:AbstractMPIReconstructionAlgorithm}, u::Array, params::LeastSquaresParameters)

  N = size(params.S, 2)
  M = div(length(params.S), N)
  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  args = toKwargs(params.solverParams)
  args[:reg] = params.reg
  args[:sparseTrafo] = params.op
  solv = createLinearSolver(params.solver, params.S; args...)

  for l=1:L
    d = solve(solv, u[:, l])
    if !isnothing(params.op)
      d[:] = params.op*d
    end
    c[:, l] = real(d)
  end

  return c

end