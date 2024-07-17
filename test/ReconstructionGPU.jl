for arrayType in arrayTypes
  @testset "Single Patch Reconstruction: $arrayType" begin
    bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
    b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))
    names = (:color, :x, :y, :z, :time)
    values = (1:1,
        -27.375u"mm":1.25u"mm":11.375u"mm",
        -11.375u"mm":1.25u"mm":27.375u"mm",
        0.0u"mm":1.0u"mm":0.0u"mm",
        0.0u"ms":0.6528u"ms":0.0u"ms")
  
    # standard reconstruction
    params = Dict{Symbol, Any}()
    params[:SNRThresh] = 5
    params[:frames] = 1:1
    params[:minFreq] = 80e3
    params[:recChannels] = 1:2
    params[:iterations] = 1
    params[:spectralLeakageCorrection] = true
    params[:sf] = bSF
    params[:reg] = [L2Regularization(0.0f0)]
    params[:solver] = CGNR
  
    c1 = reconstruct("SinglePatch", b; params...)
  
    c2 = reconstruct("SinglePatch", b; params..., arrayType = arrayType)

    @test isapprox(arraydata(c1), arraydata(c2))
  end

  @testset "Sparse Single Patch Reconstruction: $arrayType" begin
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
    params[:solver] = CGNR
    params[:redFactor] = redFactor
  
    c1 = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "FFT", useDFFoV = false)
  
    c2 = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "FFT", useDFFoV = false, arrayType = arrayType)

    @test isapprox(arraydata(c1), arraydata(c2))
  end
end