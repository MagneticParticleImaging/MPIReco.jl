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

    db = DMPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"), worker = 1)
    c2d = reconstruct("SinglePatch", db; params..., arrayType = arrayType)
    @test isapprox(arraydata(c2), arraydata(c2d))

    params[:weightingParams] = WhiteningWeightingParameters(whiteningMeas = bSF)
    c3 = reconstruct("SinglePatch", b; params...)
    c4 = reconstruct("SinglePatch", b; params..., arrayType = arrayType)
    @test isapprox(arraydata(c3), arraydata(c4))
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
    params[:sparseTrafo] = "FFT"
    params[:useDFFoV] = false
  
    c1 = reconstruct("SinglePatchSparse", b; params...)
    c2 = reconstruct("SinglePatchSparse", b; params..., arrayType = arrayType)
    @test isapprox(arraydata(c1), arraydata(c2))

    params[:weightingParams] = WhiteningWeightingParameters(whiteningMeas = bSF)
    c3 = reconstruct("SinglePatchSparse", b; params...)
    c4 = reconstruct("SinglePatchSparse", b; params..., arrayType = arrayType)
    @test isapprox(arraydata(c3), arraydata(c4))
  end

  arrayType == JLArray || @testset "Multi Patch Reconstruction: $arrayType" begin
    dirs = ["1.mdf", "2.mdf", "3.mdf", "4.mdf"]
    b = MultiMPIFile(joinpath.(datadir, "measurements", "20211226_203916_MultiPatch", dirs))
  
    dirs = ["8.mdf", "9.mdf", "10.mdf", "11.mdf"]
    bSFs = MultiMPIFile(joinpath.(datadir, "calibrations", dirs))

  
    mapping = [1,2,3,4]
    freq = filterFrequencies(bSFs, SNRThresh=5, minFreq=80e3)
    S = [getSF(SF,freq,nothing,Kaczmarz, bgcorrection=false)[1] for SF in bSFs]
    SFGridCenter = zeros(3,4)
  
    FFPos = zeros(3,4)
    FFPos[:,1] = [-0.008, 0.008, 0.0]
    FFPos[:,2] = [-0.008, -0.008, 0.0]
    FFPos[:,3] = [0.008, 0.008, 0.0]
    FFPos[:,4] = [0.008, -0.008, 0.0]
    ffPos = CustomFocusFieldPositions(FFPos)
  

    params = Dict{Symbol, Any}()
    params[:SNRThresh] = 5
    params[:frames] = [1]
    params[:minFreq] = 80e3
    params[:recChannels] = 1:2
    params[:iterations] = 50
    params[:spectralLeakageCorrection] = false
    params[:sf] = bSFs
    params[:reg] = [L2Regularization(0.0f0)]
    params[:tfCorrection] = false
    params[:solver] = CGNR

    opParams = ExplicitMultiPatchParameter(;tfCorrection = false, systemMatrices = S, SFGridCenter = SFGridCenter, mapping = mapping)
    params[:opParams] = opParams
    params[:ffPos] = ffPos
    params[:ffPosSF] = ffPos

    c1 = reconstruct("MultiPatch", b; params...)
    c2 = reconstruct("MultiPatch", b; params..., arrayType = arrayType)
    @test isapprox(arraydata(c1), arraydata(c2), rtol = 0.02)

    c3 = reconstruct("MultiPatch", b; params..., weightingParams = RowNormWeightingParameters())
    c4 = reconstruct("MultiPatch", b; params..., weightingParams = RowNormWeightingParameters(), arrayType = arrayType)
    @test isapprox(arraydata(c3), arraydata(c4), rtol = 0.02)

    # Test Kaczmarz since it uses specific functions and not just mul!
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
    c5 = reconstruct("MultiPatch", b; params...)
    c6 = reconstruct("MultiPatch", b; params..., arrayType = arrayType)
    @test isapprox(arraydata(c5), arraydata(c6))
  end

end