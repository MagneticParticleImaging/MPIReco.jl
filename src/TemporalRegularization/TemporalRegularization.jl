

mutable struct TemporalRegularizationOperator{V<:AbstractMatrix}
  S::V
  Φ::V
  N::Int
  M::Int
  L::Int # Number of all raw data frames
  idxFG::Vector{Int} # Sampling indices of foreground (sorted)
  idxBG::Vector{Int} # Sampling indices of background (sorted)
  idxCoeffsFG::Vector{Vector{Int}}
  idxCoeffsBG::Vector{Vector{Int}}
  coeffsFG::Vector{Vector{Float32}}
  coeffsBG::Vector{Vector{Float32}}
end



function reconstructionTempReg(S, bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array, grid;
  frames = nothing, bEmpty = nothing, emptyMeas= bEmpty, bgFrames = 1, numAverages = 1,
  bgDict = nothing, spectralLeakageCorrection=false, λ = 0.001, β = 0.01, solver = "kaczmarz", 
  bgCorrectionInternal=false, interpMeth=:linear, idxFG=nothing, idxBG=nothing, kargs...) where {T<:MPIFile}

  bgCorrection = emptyMeas != nothing ? true : false

  @debug "Loading emptymeas ..."
  if emptyMeas != nothing
      uEmpty = getMeasurementsFD(emptyMeas, frequencies=freq, frames=bgFrames, numAverages=length(bgFrames),
       spectralLeakageCorrection=spectralLeakageCorrection, bgCorrection=bgCorrectionInternal)
  end

  frames == nothing && (frames = 1:acqNumFrames(bMeas))

  u = getMeasurementsFD(bMeas, frequencies=freq, frames=frames, numAverages=numAverages, spectralLeakageCorrection=spectralLeakageCorrection,
                          bgCorrection=bgCorrectionInternal)

  L = -fld(-length(frames),numAverages) # number of tomograms to be reconstructed

  idxFG == nothing && (idxFG = collect(1:L))
  idxBG == nothing && (idxBG = collect(1:L))

  if emptyMeas != nothing
    u = u .- uEmpty
  end

  @info size(u)  size(S) size(bgDict) L

  #initialize output
  #image = initImage(bSF,bMeas,L,numAverages,grid,false)

  if sum(abs.(λ)) > 0 
    trace = calculateTraceOfNormalMatrix(S,nothing)
    λ *= trace / size(S,1)
  end

  Op = TemporalRegularizationOperator(Float32(1/(sqrt(λ)))*S, Float32(1/(sqrt(β)))*bgDict, L; 
                                      interpMeth, idxFG, idxBG )

  global Top = Op

  c = reconstructionTempReg(Op, u; λ, β, solver, shape=shape(grid), kargs...)

  cInterp = similar(c, size(c,1), Op.L)
  for l=1:Op.L
    for κ=1:length(Op.idxCoeffsFG[l])
      cInterp[:,l] .+= Op.coeffsFG[l][κ] * c[:, Op.idxCoeffsFG[l][κ]]
    end
  end

  image = initImage(bSF, bMeas, Op.L, numAverages, grid, false)

  @info size(image)

  image[:] .= vec(cInterp)

  @info size(image)

  return image 
end

global Top = nothing


"""
Joint reconstruction of MPI signal and background. Implemented in a low-level
fashion
"""
function reconstructionTempReg(Op::TemporalRegularizationOperator, u::Array;
                        β = 0.1, λ=0.001, solver = "kaczmarz", enforceReal=false, enforcePositive=false,
                        relativeLambda=true, backgroundCoefficients = nothing, kargs...)

  MSub, NSub, Q, J, Γ = getSizes(Op)

  u = reshape(u, :, Op.L)
  c = zeros(NSub, J)

  constraintMask = zeros(Bool, Op.N)
  constraintMask[1:(NSub*J)] .= 1

  solv = createLinearSolver(solver, Op; λ=1.0, constraintMask=constraintMask, 
                            enforceReal=enforceReal,
			                      enforcePositive=enforcePositive, kargs...)

  @info "SIZES" size(Op) size(u) size(reshape(u,:,1))


  y = solve(solv, vec(u))
  c = real.( reshape(y[1:(NSub*J),:],NSub,J) ./ sqrt(λ) ) 

  #if backgroundCoefficients != nothing
  #  append!(backgroundCoefficients, vec(y[(N+1):end,:] ./ sqrt(β)))
  #end

  return c
end





getSizes(Op::TemporalRegularizationOperator) = 
      Op.M÷Op.L, size(Op.S,2), size(Op.Φ,2), length(Op.idxFG), length(Op.idxBG) 
 
Base.eltype(TempOp::TemporalRegularizationOperator) = eltype(TempOp.S)

function size(TempOp::TemporalRegularizationOperator, i::Int)
  if i==2
    return TempOp.N
  elseif i==1
    return TempOp.M
  else
    error("bounds error")
  end
end
size(TempOp::TemporalRegularizationOperator) = (TempOp.M,TempOp.N)
length(TempOp::TemporalRegularizationOperator) = TempOp.M*TempOp.N

function TemporalRegularizationOperator(S, Φ, L; interpMeth=:linear, idxFG=1:L, idxBG=1:L )

  # check if first and last frame are sampled
  if idxFG[1] != 1 || idxFG[end] != L ||
    idxBG[1] != 1 || idxBG[end] != L
    error("idxFG or idxBG do not contain first and last frame!")
  end

  idxCoeffsFG, idxCoeffsBG, coeffsFG, coeffsBG = 
       calculateInterpolatingCoefficients(L, idxFG, idxBG, interpMeth)

  N = size(S,2)*length(idxFG)+size(S,2)*length(idxBG)
  M = size(S,1)*L

  @info "gluck" N M size(S)

  return TemporalRegularizationOperator(S, Φ, N, M, L, idxFG, idxBG, 
                           idxCoeffsFG, idxCoeffsBG, coeffsFG, coeffsBG)
end


function calculateInterpolatingCoefficients(L, idxFG, idxBG, interpMeth)
  idxCoeffsFG = Vector{Vector{Int}}(undef,0)
  idxCoeffsBG = Vector{Vector{Int}}(undef,0)
  coeffsFG = Vector{Vector{Float32}}(undef,0)
  coeffsBG = Vector{Vector{Float32}}(undef,0)

  if interpMeth==:linear
    function help(idx, idxCoeffs, coeffs)
      l = 1
      for j=1:length(idx)
        r = idx[j]:idx[min(j+1,end)]
        #push!(idxCoeffs, [idx[j]])
        push!(idxCoeffs, [j])
        push!(coeffs, [1.0])
        l += 1
        for q=2:length(r)-1
          α = (q-1)/(length(r)-1)
          #push!(idxCoeffs, [idx[j],idx[j+1]] )
          push!(idxCoeffs, [j,j+1] )
          push!(coeffs, [1-α,α] )
          l += 1
        end
      end
    end
    help(idxFG, idxCoeffsFG, coeffsFG)
    help(idxBG, idxCoeffsBG, coeffsBG)
  else
    error("Interpolation method $(interpMeth) is right now not supported")
  end
  return idxCoeffsFG, idxCoeffsBG, coeffsFG, coeffsBG
end


### The following is intended to use the standard kaczmarz method ###

function calculateTraceOfNormalMatrix(Op::TemporalRegularizationOperator, weights)
  return calculateTraceOfNormalMatrix(Op.S, weights)
end

setlambda(::TemporalRegularizationOperator, ::Any) = nothing

getSubIdxRows(k,M) = mod1(k,M), (k-1)÷M+1

function RegularizedLeastSquares.dot_with_matrix_row(Op::TemporalRegularizationOperator, 
                                                     x::AbstractArray{T}, k::Integer) where T
  
  MSub, NSub, Q, J, Γ = getSizes(Op)

  #@info "Blub"
  #@info MSub, NSub, Q, J, Γ 
  #@info size(x), k

  i,l = getSubIdxRows(k, MSub)
  xFG = reshape( view(x, 1:NSub*J), NSub, J)
  xBG = reshape( view(x, (NSub*J+1):(NSub*J+Q*Γ)), Q, Γ)
  
  #@info l 
  #@info size(Op.idxCoeffsFG)
  #@info Op.idxCoeffsFG[l]

  tmpFG = Op.coeffsFG[l][1] * xFG[:, Op.idxCoeffsFG[l][1]]
  for κ=2:length(Op.idxCoeffsFG[l])
    tmpFG .+= Op.coeffsFG[l][κ] * xFG[:, Op.idxCoeffsFG[l][κ]]
  end

  tmpBG = Op.coeffsBG[l][1] * xBG[:, Op.idxCoeffsBG[l][1]]
  for κ=2:length(Op.idxCoeffsBG[l])
    tmpBG .+= Op.coeffsBG[l][κ] * xBG[:, Op.idxCoeffsBG[l][κ]]
  end

  tmp = zero(T)
  @simd  for n = 1:NSub
    @inbounds tmp += Op.S[i,n]*tmpFG[n] 
  end
  @simd  for q = 1:Q
    @inbounds tmp += Op.Φ[i,q]*tmpBG[q]  
  end
  return tmp
end


function RegularizedLeastSquares.kaczmarz_update!(Op::TemporalRegularizationOperator, 
                                                  x::AbstractArray, k::Integer, β)
  MSub, NSub, Q, J, Γ = getSizes(Op)
  i, l = getSubIdxRows(k, MSub)
  xFG = reshape( view(x, 1:NSub*J), NSub, J)
  xBG = reshape( view(x, (NSub*J+1):(NSub*J+Q*Γ)), Q, Γ)
  
  for κ=1:length(Op.idxCoeffsFG[l])
    α = β * conj(Op.coeffsFG[l][κ])
    #@simd 
    for n = 1:NSub
      #@inbounds 
      xFG[n,Op.idxCoeffsFG[l][κ]] += α * conj(Op.S[i,n])
    end
  end


  for κ=1:length(Op.idxCoeffsBG[l])
    α = β * conj(Op.coeffsBG[l][κ])
    #@simd 
    for q = 1:Q
     # @info q Q
      #@inbounds 
      xBG[q,Op.idxCoeffsBG[l][κ]] += α * conj(Op.Φ[i,q])
    end
  end
  return
end

function RegularizedLeastSquares.initkaczmarz(Op::TemporalRegularizationOperator,λ,weights::Vector)
  T = typeof(real(Op.S[1]))
  denom = T[]
  rowindex = Int64[]
  MSub, NSub, Q, S, Γ = getSizes(Op)

  s² = [rownorm²(Op.S,i) for i=1:MSub]
  φ² = [rownorm²(Op.Φ,i) for i=1:MSub]

  for l=1:Op.L
    for i=1:MSub
      if s²[i] > 0 || φ²[i] > 0
        k = i+MSub*(l-1)
        push!(denom,1 / (s²[i]*sum(abs.(Op.coeffsFG[l]).^2) + 
                         φ²[i]*sum(abs.(Op.coeffsBG[l]).^2) + λ)) 
        push!(rowindex,k) 
      end
    end
  end

  denom, rowindex
end