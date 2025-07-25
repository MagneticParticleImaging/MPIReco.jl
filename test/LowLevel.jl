@testset "LowLevel" begin
  
  @testset "Single-Patch" begin
    bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
    b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))
    params = Dict{Symbol, Any}()
    params[:SNRThresh] = 5
    params[:frames] = 1:1
    params[:minFreq] = 80e3
    params[:recChannels] = 1:2
    params[:iterations] = 1
    params[:spectralLeakageCorrection] = true
    params[:sf] = bSF
    params[:reg] = [L2Regularization(0.0f0)]
    params[:solver] = Kaczmarz

    # No weighting
    high = build(MPIRecoPlan("SinglePatch"; params...))
    cHigh = reconstruct(high, b)
    
    S = high.S
    u = process(high, high.params.pre, b, high.freqs)
    cLow = reconstruct("LowLevel", u; S = S, iterations = params[:iterations], reg = params[:reg], solver = params[:solver])
    
    cLow = reshape(cLow, size(cHigh))
    @test isapprox(cHigh.data.data, cLow)

    # Weighted
    params[:weightingParams] = WhiteningWeightingParameters(whiteningMeas = bSF)
    high = build(MPIRecoPlan("SinglePatch"; params...))
    cHigh = reconstruct(high, b)
    
    S = high.S
    u = process(high, high.params.pre, b, high.freqs)
    weights = high.weights
    cLow = reconstruct("LowLevel", u; S = S, iterations = params[:iterations], reg = params[:reg], solver = params[:solver], weights = weights)
    
    cLow = reshape(cLow, size(cHigh))
    @test isapprox(cHigh.data.data, cLow)
  end

    @testset "Single-Patch Sparse" begin
    bSF = MPIFile(joinpath(datadir, "calibrations", "7.mdf"))
    b = MPIFile(joinpath(datadir, "measurements","20211226_204612_Dice", "1.mdf"))
    redFactor = 0.01
    params = Dict{Symbol, Any}()
    params[:SNRThresh] = 2
    params[:frames] = 1:100
    params[:minFreq] = 80e3
    params[:recChannels] = 1:2
    params[:iterations] = 3
    params[:spectralLeakageCorrection] = false
    params[:sf] = bSF
    params[:reg] = [L2Regularization(0.1f0)]
    params[:numAverages] = 100
    params[:solver] = Kaczmarz
    params[:redFactor] = redFactor

    # No weighting
    high = build(MPIRecoPlan("SinglePatchSparse"; params..., sparseTrafo = "FFT", useDFFoV = false))
    cHigh = reconstruct(high, b)
    
    S = high.S
    u = process(high, high.params.pre, b, high.freqs)
    op = MPIReco.getLinearOperator(high, high.params.reco)
    cLow = reconstruct("LowLevel", u; S = S, op = op, iterations = params[:iterations], reg = params[:reg], solver = params[:solver])
    
    cLow = reshape(cLow, size(cHigh))
    @test isapprox(cHigh.data.data, cLow)

    # Weighted
    params[:weightingParams] = WhiteningWeightingParameters(whiteningMeas = bSF)
    high = build(MPIRecoPlan("SinglePatchSparse"; params..., sparseTrafo = "FFT", useDFFoV = false))
    cHigh = reconstruct(high, b)
    
    S = high.S
    u = process(high, high.params.pre, b, high.freqs)
    op = MPIReco.getLinearOperator(high, high.params.reco)
    weights = high.weights
    cLow = reconstruct("LowLevel", u; S = S, op = op, iterations = params[:iterations], reg = params[:reg], solver = params[:solver], weights = weights)
    
    cLow = reshape(cLow, size(cHigh))
    @test isapprox(cHigh.data.data, cLow)
  end

  @testset "Multi-Patch" begin
    bSF = MultiMPIFile([joinpath(datadir, "calibrations", "12.mdf")])
    dirs = ["1.mdf", "2.mdf", "3.mdf", "4.mdf"]
    b = MultiMPIFile(joinpath.(datadir, "measurements", "20211226_203916_MultiPatch", dirs))
    params = Dict{Symbol, Any}()
    params[:SNRThresh] = 5
    params[:frames] = [1]
    params[:minFreq] = 80e3
    params[:recChannels] = 1:2
    params[:iterations] = 1
    params[:spectralLeakageCorrection] = false
    params[:sf] = bSF
    params[:reg] = [L2Regularization(0.0f0)]
    params[:tfCorrection] = false
    params[:solver] = Kaczmarz

    # No weighting
    high = build(MPIRecoPlan("MultiPatch"; params...))
    cHigh = reconstruct(high, b)
    
    S = copy(high.ffOp)
    u = process(high, high.params.pre, b, high.freqs)
    cLow = reconstruct("LowLevel", u; S = S, iterations = params[:iterations], reg = params[:reg], solver = params[:solver])
    
    cLow = reshape(cLow, size(cHigh))
    @test isapprox(cHigh.data.data, cLow)

    # Weighted
    params[:weightingParams] = RowNormWeightingParameters()
    high = build(MPIRecoPlan("MultiPatch"; params...))
    cHigh = reconstruct(high, b)
    
    S = copy(high.ffOp)
    u = process(high, high.params.pre, b, high.freqs)
    weights = high.weights
    cLow = reconstruct("LowLevel", u; S = S, iterations = params[:iterations], reg = params[:reg], solver = params[:solver], weights = weights)
    
    cLow = reshape(cLow, size(cHigh))
    @test isapprox(cHigh.data.data, cLow)
  end

end