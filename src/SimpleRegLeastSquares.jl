# TODO Maybe chose better names depending on name of other least squares implementation, (for example hands-free)
mutable struct SimpleRegularizedLeastSquares <: AbstractMPIReconstructionAlgorithm
  # Parameter
  S::Any
  sparseTrafo = nothing # TODO Type
  solver::AbstractLinearSolver
  # "Algorithm" fields
  outputChannel::Channel{Any}
end

function SimpleRegularizedLeastSquares(S; sparseTrafo=nothing,
  lambd=0.0, lambda=lambd, λ=lambda, solver="kaczmarz",
  weights=nothing, enforceReal=true, enforcePositive=true,
  relativeLambda=true, kargs...)

  N = size(S, 2)

  if sum(abs.(λ)) > 0 && solver != "fusedlasso" && relativeLambda
    trace = calculateTraceOfNormalMatrix(S, weights)
    if isa(λ, AbstractVector)
      λ[1:1] *= trace / N
    else
      λ *= trace / N
    end
    setlambda(S, λ)
  end

  solv = createLinearSolver(solver, S; weights=weights, λ=λ,
    sparseTrafo=sparseTrafo, enforceReal=enforceReal,
    enforcePositive=enforcePositive, kargs...)
  return SimpleRegularizedLeastSquares(S, sparseTrafo, solv, Channel{Any}(32))
end

function put!(algo::SimpleRegularizedLeastSquares, u::Array)
  S = algo.S
  solv = algo.solver
  sparseTrafo = algo.sparseTrafo

  N = size(S, 2) #prod(shape)
  M = div(length(S), N)

  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  for l = 1:L

    d = solve(solv, u[:, l])

    if sparseTrafo != nothing
      d[:] = sparseTrafo * d #backtrafo from dual space
    end

    #if typeof(B)==LinearSolver.DSTOperator
    #	d=onGridReverse(d,shape)
    #end
    c[:, l] = real(d) # this one is allocating
  end
  put!(algo.outputChannel, c)
end

take!(algo::SimpleRegularizedLeastSquares) = take!(algo.outputChannel)