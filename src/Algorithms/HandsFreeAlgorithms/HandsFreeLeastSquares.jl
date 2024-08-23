export HandsFreeSolverParameters
Base.@kwdef struct HandsFreeSolverParameters <: AbstractSolverParameters{Kaczmarz}
  iterbounds::Tuple{Int64, Int64}=(1, 25)
  enforceReal::Bool=true
  enforcePositive::Bool=true
  normalizeReg::AbstractRegularizationNormalization = SystemMatrixBasedNormalization()
  startλ::Float64=5.0
  SNRbounds::Tuple{Float64, Float64}=(0, 1)
  stoppingParas::Tuple{Float64, Float64}=(0.25, 2.0)
  equalizeIters::Bool=false
  flattenIters::Bool=false
end

function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::LeastSquaresParameters{Kaczmarz, O, SF, R, SL, W}, u::AbstractArray, snr::AbstractVector) where {O, SF, R, SL <: HandsFreeSolverParameters, W}

  solverParams = params.solverParams
  N = size(params.S, 2)
  M = div(length(params.S), N)
  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  reg, _ = prepareRegularization([L2Regularization(real(eltype(u))(solverParams.startλ))], params)

  S = params.S
  if !isnothing(params.weights)
    S = ProdOp(WeightingOp(params.weights), S)
    for l = 1:L
      u[:, l] = params.weights.*u[:, l]
    end
  end

  it = zeros(L)
  w_it = zeros(L)
  curv = [zeros(solverParams.iterbounds[2]) for i=1:L]

  solv = Kaczmarz(S, snr; reg = reg, iterbounds = solverParams.iterbounds, stoppingParas = solverParams.stoppingParas, SNRbounds = solverParams.SNRbounds, normalizeReg = solverParams.normalizeReg)

  for l=1:L
    d = solve!(solv, u[:, l])

    # Update curve meta data
    curv[l] = solv.state.curvatures
    it[l] = solv.state.iteration
    w_it[l] = solv.state.wanted_iters

    # Store frame
    if !isnothing(params.op)
      d[:] = params.op*d
    end
    c[:, l] = Array(real(d))
    
    if solverParams.equalizeIters && l >= 3
      solv.state.expected_iters = round(Int, sum(it[l-2:l].*[0.1, 0.3, 0.6]))
    end
  end

  if solverParams.flattenIters
    for l=1:L
      iter = round(Int, mean(it[maximum((1, l-5)):minimum((l+5, L))]))
      solv.state.iterbounds = (iter, iter)
      solv.state.expected_iters = iter

      d = solve!(solv, u[:, l])

      # Update curve meta data
      curv[l] = solv.state.curvatures
      it[l] = solv.state.iteration
      w_it[l] = solv.state.iteration
  
      # Store frame
      if !isnothing(params.op)
        d[:] = params.op*d
      end
      c[:, l] = Array(real(d))  
    end
  end

  return c

end