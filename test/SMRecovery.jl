using MPIReco

@testset "system matrix recovery" begin
  # patterns of system matrix
  Random.seed!(1234)
  bSF = MPIFile("$datadir/SF_MP")
  shape = calibSize(bSF)
  numFreqs = 3
  f = filterFrequencies(bSF, sortBySNR=true)
  f = f[shuffle(1:100)[1:numFreqs]]
  S = getSystemMatrix(bSF,bgCorrection=true)
  S = S[:,f]

  # undersampled measurement
  numSamples = Int64(div(prod(shape),2))
  samplingIdx = shuffle(1:prod(shape))[1:numSamples]
  y = S[samplingIdx,:]

  # solver parameters
  params = Dict{Symbol,Any}()
  params[:shape] = Tuple(shape)
  params[:ρ_l1] = 2.0
  params[:λ_l1] = 2.e-2
  params[:iterationsInner] = 50
  params[:iterations] = 10
  params[:relTol] = 1.e-2
  params[:absTol] = 1.e-4

  # DCT based recovery
  @info "recovery using DCT"
  S_dct = smRecovery(y,samplingIdx,params)
  for k=1:numFreqs
    @test norm(S[:,k] - S_dct[:,k]) / norm(S[:,k]) ≈ 0 atol=0.11
  end

  # DCT and low rank
  @info "recovery using DCT and LR"
  params[:λ_lr] = 2.e-2
  params[:ρ_lr] = 2.e-1
  S_dctlr = smRecovery(y,samplingIdx,params)
  for k=1:numFreqs
    @test norm(S[:,k] - S_dct[:,k]) / norm(S[:,k]) ≈ 0 atol=0.11
  end
end
