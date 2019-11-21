
function getBackgroundDictionary(fSF::MPIFile, frequencies::Vector, bgDictSize::Int=2)
  idxBGFrames = measBGFrameIdx(fSF)

  D = measData(fSF, idxBGFrames)
  D_ = reshape(D, size(D,1), size(D,2)*size(D,3), size(D,4))
  bgdata = reshape(D_[:,frequencies,:],length(idxBGFrames),:)

  U,S,V = svd(transpose(bgdata))

  return U[:,1:bgDictSize]
end

getBackgroundDictionary(::Any, ::Vector, bgDictSize::Nothing) = nothing



"""
Joint reconstruction of MPI signal and background. Implemented in a low-level
fashion
"""
function reconstruction(S, u::Array, bgDict::AbstractMatrix;
                        sparseTrafo = nothing, bgDictSize = 2, beta = 0.1, β=beta,
                        lambd=0.0, lambda=lambd, λ=lambda, progress=nothing,
                        solver = "kaczmarz",
                        weights=nothing, enforceReal=false, enforcePositive=false,
                        relativeLambda=true, kargs...)

  N = size(S,2)
  M = div(length(S), N)

  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N,L)

  β /= λ

  if sum(abs.(λ)) > 0 && solver != "fusedlasso" && relativeLambda
    trace = calculateTraceOfNormalMatrix(S,weights)
    λ *= trace / N
    setlambda(S,λ)
  end

  G = transpose( cat(transpose(S), Float32(1/β)*transpose(bgDict), dims=1) )

  constraintMask = zeros(Bool, size(G,2))
  constraintMask[1:N] .= 1

  solv = createLinearSolver(solver, G; weights=weights, λ=λ, constraintMask=constraintMask,
                            sparseTrafo=sparseTrafo, enforceReal=enforceReal,
			                enforcePositive=enforcePositive, kargs...)

  progress==nothing ? p = Progress(L, 1, "Reconstructing data...") : p = progress
  for l=1:L

    d = solve(solv, u[:,l])[1:N,:]

    if sparseTrafo != nothing
      d[:] = sparseTrafo*d # backtrafo from dual space
    end

    c[:,l] = real( d )
    next!(p)
    sleep(0.001)
  end

  return c
end
