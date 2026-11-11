using MPIReco: process, AbstractMPIRecoAlgorithm

@testset "Weighting" begin
  # Weights are computed with the signature used by the reconstruction algorithms:
  # process(algoT, weightingParams, freqs, sf, S, u, arrayType)
  algoT = AbstractMPIRecoAlgorithm
  S = ComplexF32[1 2 0; 0 3 4; 5 0 0; 1 1 1] # 4 frequency components x 3 voxels
  freqs = [CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(1, 2), CartesianIndex(2, 2)]
  rowNormWeights = [1 / sqrt(sum(abs2, S[r, :])) for r in 1:size(S, 1)]

  @testset "NoWeighting" begin
    @test isnothing(process(algoT, NoWeightingParameters(), freqs, nothing, S, nothing, Array))
  end

  @testset "ChannelWeighting" begin
    w = process(algoT, ChannelWeightingParameters([0.5, 2.0]), freqs, nothing, S, nothing, Array)
    @test w ≈ [0.5, 0.5, 2.0, 2.0]
  end

  @testset "RowNormWeighting" begin
    w = process(algoT, RowNormWeightingParameters(), freqs, nothing, S, nothing, Array)
    @test w ≈ rowNormWeights
  end

  @testset "CompositeWeighting" begin
    params = CompositeWeightingParameters([ChannelWeightingParameters([0.5, 2.0]), RowNormWeightingParameters()])
    w = process(algoT, params, freqs, nothing, S, nothing, Array)
    @test w ≈ [0.5, 0.5, 2.0, 2.0] .* rowNormWeights
  end

  @testset "SMWhiteningWeighting" begin
    # Synthetic system matrix in the processed calibration layout with background frames
    N, K, C = 12, 9, 2
    isBG = fill(false, N)
    isBG[[1, 2, 6, 12]] .= true
    data = ComplexF32[sin(n + 2k + 3c) + im * cos(5n - k + c) for n = 1:N, k = 1:K, c = 1:C, j = 1:1]
    sf = defaultMDFv2InMemory()
    rxNumChannels(sf, C)
    rxNumSamplingPoints(sf, 2(K - 1))
    rxBandwidth(sf, 1.0e6)
    acqNumFrames(sf, N)
    acqNumPeriodsPerFrame(sf, 1)
    measDataRaw(sf, data)
    measIsFourierTransformed(sf, true)
    measIsFastFrameAxis(sf, true)
    measIsFramePermutation(sf, true)
    measIsFrequencySelection(sf, false)
    measIsBGFrame(sf, isBG)
    measIsBGCorrected(sf, false)
    measIsTFCorrected(sf, false)
    measIsSpectralLeakageCorrected(sf, false)
    measIsSparsityTransformed(sf, false)

    noise = std(data[isBG, :, :, 1], dims = 1)[1, :, :][freqs]
    w = process(algoT, SMWhiteningWeightingParameters(), freqs, sf, S, nothing, Array)
    @test w ≈ minimum(noise) ./ noise
    @test maximum(w) ≈ 1

    tf = ComplexF64[1 + k / K + im * c / C for k = 1:K, c = 1:C]
    rxTransferFunction(sf, tf)
    noiseTF = noise ./ abs.(tf[freqs])
    # tfCorrection defaults to rxHasTransferFunction(sf)
    w = process(algoT, SMWhiteningWeightingParameters(), freqs, sf, S, nothing, Array)
    @test w ≈ minimum(noiseTF) ./ noiseTF
    w = process(algoT, SMWhiteningWeightingParameters(tfCorrection = false), freqs, sf, S, nothing, Array)
    @test w ≈ minimum(noise) ./ noise
  end
end
