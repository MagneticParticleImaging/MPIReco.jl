

function reconstructionTempRegTikhonov(S, bSF::Union{T,Vector{T}}, bMeas::MPIFile, freq::Array, grid;
  frames = nothing, bEmpty = nothing, emptyMeas= bEmpty, bgFrames = 1, numAverages = 1,
  bgDict = nothing, spectralLeakageCorrection=false, λ = [0.001,0.001], β = [0.01,0.01], solver = "kaczmarz", 
  bgCorrectionInternal=false, kargs...) where {T<:MPIFile}

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

  if emptyMeas != nothing
    u = u .- uEmpty
  end

  @info size(u)  size(S) size(bgDict) L

  if sum(abs.(λ)) > 0 
    trace = calculateTraceOfNormalMatrix(S,nothing)
    λ *= trace / size(S,2)
  end

  c = reconstructionTempRegTikhonov(S, bgDict, u; λ, β, solver, shape=shape(grid), kargs...)

  image = initImage(bSF, bMeas, L, numAverages, grid, false)

  @info size(image)

  image[:] .= vec(c)

  @info size(image)

  return image 
end

"""
Joint reconstruction of MPI signal and background. Implemented in a low-level
fashion
"""
function reconstructionTempRegTikhonov(S, Φ, u::Array;
                        λ = [0.001,0.001], β = [0.01,0.01], solver = "kaczmarz", 
                        enforceReal=false, enforcePositive=false,
                        relativeLambda=true, backgroundCoefficients = nothing, kargs...)

  N = size(S, 2)
  M = div(length(S), N)
  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  λf = Float32.(1.0 ./ sqrt.(λ))
  βf = Float32.(1.0 ./ sqrt.(β))

  G1 = [λf[1]*S zeros(Float32,size(S)); λf[2]*S  λf[2]*S]
  G2 = [βf[1]*Φ zeros(Float32,size(Φ)); βf[2]*Φ  βf[2]*Φ]
  #G1 = [zeros(Float32,size(S)) λf[1]*S; λf[2]*S  λf[2]*S]
  #G2 = [zeros(Float32,size(Φ)) βf[1]*Φ; βf[2]*Φ  βf[2]*Φ]
  G = transpose( cat(transpose(G1), transpose(G2), dims=1) )

  constraintMask = zeros(Bool, size(G, 2))
  constraintMask[1:N] .= 1

  solv = createLinearSolver(solver, G; λ=1.0, constraintMask=constraintMask, 
                            enforceReal=enforceReal,
			                      enforcePositive=enforcePositive, kargs...)

  for l=1:L-1
    w = cat(u[:,l],u[:,l+1],dims=1)
    y = solve(solv, w)
    #c[:,l] = real.( y[N+1:2*N] .* λf[1]) 
    c[:,l] = real.( y[1:N] .* λf[1] + 0.5*(y[N+1:2*N] .* λf[2])) 

  end
  #if backgroundCoefficients != nothing
  #  append!(backgroundCoefficients, vec(y[(N+1):end,:] ./ sqrt(β)))
  #end

  return c
end

