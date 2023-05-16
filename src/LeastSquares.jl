export LeastSquaresParameters

abstract type AbstractRegularizationParameters <: AbstractMPIRecoParameters end
abstract type AbstractSolverIterationParameters <: AbstractMPIRecoParameters end

export LeastSquaresParameters
Base.@kwdef struct LeastSquaresParameters{L<:AbstractLinearSolver, O, M, R<:AbstractRegularizationParameters, P<:AbstractSolverIterationParameters} <: AbstractMPIRecoParameters
  solver::Type{L}
  op::O
  S::M
  reg::R = L2Regularization(0.1)
  solvP::P
end

# TODO place weights and more

export SimpleSolverIterationParameters
Base.@kwdef struct SimpleSolverIterationParameters <: AbstractSolverIterationParameters
  iterations::Int64=10
  enforceReal::Bool=false
  enforcePositive::Bool=true
end

export SimpleRegularizationParameters
Base.@kwdef struct SimpleRegularizationParameters <: AbstractRegularizationParameters
  λ::Vector{Float64} = [0.1]
  regName::Vector{String} = ["L2"]
end
export L2Regularization
L2Regularization(λ::Float64) = SimpleRegularizationParameters([λ], ["L2"])

function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, u::Array, params::LeastSquaresParameters)

  N = size(params.S, 2)
  M = div(length(params.S), N)
  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  # TODO Where to solve this? Ideally here I just want to have the creation and interaction with the solver, everything else is already done
  #if sum(abs.(λ)) > 0 && params.solver != FusedLasso && params.relativeLambda
  #  trace = calculateTraceOfNormalMatrix(params.S,weights)
  #  if isa(λ,AbstractVector) 
  #    λ[1:1] *= trace / N
  #  else
  #    λ *= trace / N
  #  end
  #  #setlambda(S,λ) dead code?
  #end

  args = toKwargs([params.reg, params.solvP])
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