export HandsFreeKaczmarzState
mutable struct HandsFreeKaczmarzState{T, Tc <: Union{T, Complex{T}}, vecTc <: AbstractArray{Tc}, vecT <: AbstractArray{T}, R <: AbstractRegularization}  <: RegularizedLeastSquares.AbstractSolverState{Kaczmarz}
  # Usual Kaczmarz
  u::vecTc
  x::vecTc
  vl::vecTc
  regs::Vector{R}
  εw::vecT
  τl::Tc
  αl::Tc
  iteration::Int64
  # HandsFree
  freqSNRs::vecT
  curvatures::vecT
  norm_cl::vecT
  resid::vecT
  iterbounds::Tuple{Int64,Int64}
  stoppingParas::Tuple{Float64,Float64}
  expected_iters::Int64
  wanted_iters::Int64
  SNRbounds::Tuple{Float64,Float64}
end

function RegularizedLeastSquares.Kaczmarz(A, freqSNRs
                ; reg = L2Regularization(real(eltype(A))(5.0)) # encodes start value for λ
                , normalizeReg::AbstractRegularizationNormalization = NoNormalization()
                , iterbounds::Tuple{Int64,Int64}=(1,50)
                , stoppingParas::Tuple{Float64,Float64}=(0.25,2.0)
                , expected_iters::Int64=0
                , wanted_iters::Int64=0
                , SNRbounds::Tuple{Float64,Float64}=(60.0,1.5)
                )

  T = real(eltype(A))

  min_iter, max_iter = iterbounds

  # Prepare regularization terms
  reg = isa(reg, AbstractVector) ? reg : [reg]
  idx = findsink(L2Regularization, reg)
  regs = AbstractRegularization[]
  if length(idx) == max_iter # User provided enough L2Regularization terms
    regs = reg[idx]
  elseif length(idx) == 1
    startλ = λ(reg[idx])
    regs = [L2Regularization(startλ / (1 + (0.2 * i - 0.2)^5)) for i in 1:max_iter]
  else
    error("HandsFreeKaczmarz requires either one or as many as the upper itterbound L2Regularization terms, found $(length(idx))")
  end
  deleteat!(reg, idx)
  regs = RegularizedLeastSquares.normalize(Kaczmarz, normalizeReg, regs, A, nothing)

  # Tikhonov matrix is only valid with NoNormalization or SystemMatrixBasedNormalization
  if λ(first(regs)) isa AbstractVector
    error("Tikhonov matrix for Handsfree Kaczmarz is not yet implemented")
  end

  indices = findsinks(AbstractProjectionRegularization, reg)
  other = AbstractRegularization[reg[i] for i in indices]
  deleteat!(reg, indices)
  if length(reg) == 1
    push!(other, reg[1])
  elseif length(reg) > 1
    error("Kaczmarz does not allow for more than one additional regularization term, found $(length(reg))")
  end
  other = identity.(other)

  # setup denom and rowindex
  denom, rowindex = inithandsfree(A)
  rowIndexCycle = collect(1:length(rowindex))

  M,N = size(A)

  u  = zeros(eltype(A),M)
  x = zeros(eltype(A),N)
  vl = zeros(eltype(A),M)
  εw = zeros(T, max_iter)
  τl = zero(eltype(A))
  αl = zero(eltype(A))

  state = HandsFreeKaczmarzState(u, x, vl, regs, εw, τl, αl, 0, freqSNRs, zeros(T, max_iter), zeros(T, max_iter), zeros(T, max_iter), iterbounds, stoppingParas, expected_iters, wanted_iters, SNRbounds)

  return Kaczmarz(A, nothing, other, denom, rowindex, rowIndexCycle,
    false, 0, T[], false,
    Int64(0), normalizeReg, max_iter, state)
end

function inithandsfree(A)
  T = real(eltype(A))
  denom = T[]
  rowindex = Int64[]

  for i = 1:size(A, 1)
    s² = rownorm²(A,i)
    if s²>0
      push!(denom,s²) # only compute rownorm2, λ is computed during iterations
      push!(rowindex,i)
    end
  end
  return denom, rowindex
end

function RegularizedLeastSquares.init!(solver::Kaczmarz, state::HandsFreeKaczmarzState{T, Tc, vecTc, vecT}, b::vecTc; x0 = 0) where {T, Tc, vecTc, vecT}
  solver.reg = RegularizedLeastSquares.normalize(solver, solver.normalizeReg, solver.reg, solver.A, b)
  state.regs = RegularizedLeastSquares.normalize(solver, solver.normalizeReg, state.regs, solver.A, b)

  #if solver.shuffleRows || solver.randomized
  #  Random.seed!(solver.seed)
  #end
  #if solver.shuffleRows
  #  shuffle!(solver.rowIndexCycle)
  #end

  # start vector
  state.x .= x0
  state.vl .= 0

  state.resid[:] .= zero(T)
  state.norm_cl[:] .= zero(T)
  state.curvatures[:] .= zero(T)


  state.u .= b
  state.εw .= sqrt.(λ.(state.regs))
  state.iteration = 1
end

function iterate(solver::Kaczmarz{matT, Nothing}, state::HandsFreeKaczmarzState{T, Tc, vecTc, vecT} = solver.state) where {matT, T, Tc, vecTc, vecT}
  if RegularizedLeastSquares.done(solver,state) return nothing end
  iteration = state.iteration

  usedIndices = filterFrequencies(solver, state)

  for (i, row) in enumerate(usedIndices) 
    RegularizedLeastSquares.iterate_row_index(solver, state, solver.A, row, i, iteration)
  end

  for r in solver.reg
    prox!(r, state.x)
  end
  
  state.norm_cl[iteration] = norm(real(state.x))
  state.resid[iteration] = norm(- sqrt(λ(state.regs[iteration])) * state.vl[usedIndices]) / length(usedIndices)
  
  if iteration == 1
      state.curvatures[iteration] = 0
  elseif iteration == 2
      dcdr = (state.norm_cl[iteration]-state.norm_cl[iteration-1]) / (state.resid[iteration]-state.resid[iteration-1])
      state.curvatures[iteration] = ( (dcdr - 0) / (state.resid[iteration]-state.resid[iteration-1]) ) / ((1 + dcdr^2)^(3/2))
  else
      dcdr_old = (state.norm_cl[iteration-1] - state.norm_cl[iteration-2]) / (state.resid[iteration-1] - state.resid[iteration-2])
      dcdr = (state.norm_cl[iteration] - state.norm_cl[iteration-1]) / (state.resid[iteration] - state.resid[iteration-1])
      state.curvatures[iteration] = ( (dcdr - dcdr_old) / (state.resid[iteration]-state.resid[iteration-1]) ) / ((1 + dcdr^2)^(3/2))
  end

  state.iteration += 1
  return state.x, state
end

function RegularizedLeastSquares.iterate_row_index(solver::Kaczmarz, state::HandsFreeKaczmarzState, A, row, index, iteration)
  state.τl = dot_with_matrix_row(A,state.x,row)
  state.αl = 1/(solver.denom[index] + state.ɛw[iteration]^2) * (state.u[row]-state.τl-state.ɛw[iteration]*state.vl[row])
  kaczmarz_update!(A,state.x,row,state.αl)
  state.vl[row] += state.αl*state.ɛw[iteration]
end


function filterFrequencies(::Kaczmarz, state::HF) where {T, Tc, HF <: HandsFreeKaczmarzState{T, Tc}}
  iteration = state.iteration
  threshl = [j > state.SNRbounds[2] ? T(j) : T(state.SNRbounds[2]) for j in (state.SNRbounds[1]) / (1 + (0.28 * iteration - 0.28)^2)][1]
  indexl = sort(findall(x -> x > threshl, state.freqSNRs)) # potentially drop sorting
  return indexl
end

function RegularizedLeastSquares.done(::Kaczmarz,state::HF) where HF <: HandsFreeKaczmarzState
  iteration = state.iteration
  if state.expected_iters != 0 && iteration >= round(Int,state.expected_iters*3/2)    
      tmp=state.expected_iters
      @info "Would like to go more than $iteration iterations, but expect $tmp iterations."
      if !(state.wanted_iters < state.expected_iters)
          state.wanted_iters = iteration+2
      end
      return true
  elseif iteration >= 2 && iteration >= state.iterbounds[1] &&  state.curvatures[iteration] > state.stoppingParas[1] * state.norm_cl[1] && abs(state.curvatures[iteration-1]*state.stoppingParas[2]) < abs(state.curvatures[iteration])
      if iteration >= round(Int,state.expected_iters*1/2)
          tmp=state.expected_iters
          @info "Stopped after $iteration iterations. Expected $tmp iterations."
          state.wanted_iters = iteration
          return true
      else
          tmp = state.expected_iters
          tmp2 = round(Int,state.expected_iters*2/5)
          state.expected_iters = tmp2
          state.wanted_iters = iteration
          @info "Would like to stop after $iteration iterations, but expect $tmp iterations. Update expected iterations for this reco to $tmp2."
          return false
      end
  elseif iteration >= state.iterbounds[2]
      tmp = state.expected_iters
      @info "Stopped at the max iter-bound of $iteration iterations. Expected $tmp iterations."
      state.wanted_iters = iteration
      return true
  else
      return false
  end
end