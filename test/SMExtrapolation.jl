using MPIReco

@testset "SystemMatrixExtrapolation" begin
    bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
    freq = sort(filterFrequencies(bSF,SNRThresh=5,minFreq=60e3,recChannels=1:2))
    SMextr1 = extrapolateSM(bSF, freq, (5,5,0); method=1)[1]
    
    SM, grid = getSF(bSF, freq, returnasmatrix=true)
    SMextr2, extrgrid = extrapolateSM(SM, grid, (5,5,0); method=1)
    SMextr3 = extrapolateSM(SM, grid, 5)[1]
    
    @test SMextr2 == transpose(SMextr1)
    @test SMextr2 == SMextr3

    SMsizeExtr = shape(grid) .+ [10;10;0]
    @test shape(extrgrid) == SMsizeExtr
    @test prod(SMsizeExtr) == size(SMextr2,1)
    
    b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))
    plan = getPlan("Single")
    setAll!(plan, :SNRThresh, 5)
    setAll!(plan, :frames, 1:1)
    setAll!(plan, :minFreq, 80e3),
    setAll!(plan, :recChannels, 1:2)
    setAll!(plan, :iterations, 1)
    setAll!(plan, :spectralLeakageCorrection, true)
    setAll!(plan, :sf, bSF)
    setAll!(plan, :reg, [L2Regularization(0.0f0)])
    setAll!(plan, :solver, Kaczmarz)
    setAll!(plan, :gridsize, calibSize(bSF))
    setAll!(plan, :fov, calibFov(bSF))  
    c1 = reconstruct(build(plan), b)

    setAll!(plan, :gridsize, [36,36,1])
    setAll!(plan, :fov, calibFov(bSF).+[0.006,0.006,0])  
    c_extr = reconstruct(build(plan), b)
    @test size(c1[1,:,:,:,1]) .+ (4,4,0) == size(c_extr[1,:,:,:,1])
    exportImage(joinpath(imgdir, "Extrapolated1.png"), Array(c_extr[1,:,:,1,1]))
    @test compareImg("Extrapolated1.png")

    setAll!(plan, :gridsize, calibSize(bSF))
    setAll!(plan, :fov, calibFov(bSF).+[0.006,0.006,0])  
    c_extr2 = reconstruct(build(plan), b)
    @test size(c1[1,:,:,:,1]) == size(c_extr2[1,:,:,:,1])
    exportImage(joinpath(imgdir, "Extrapolated2.png"), Array(c_extr2[1,:,:,1,1]))
    @test compareImg("Extrapolated2.png")

    SFdirs = ["8.mdf", "9.mdf", "10.mdf", "11.mdf"]
    bSFs = MultiMPIFile(joinpath.(datadir, "calibrations", SFdirs))
    bdirs = ["1.mdf", "2.mdf", "3.mdf", "4.mdf"]
    b = MultiMPIFile(joinpath.(datadir, "measurements", "20211226_203916_MultiPatch", bdirs))

    plan = getPlan("MultiPatch")
    setAll!(plan, :SNRThresh, 5)
    setAll!(plan, :frames, 1:1)
    setAll!(plan, :minFreq, 80e3)
    setAll!(plan, :recChannels, 1:2)
    setAll!(plan, :iterations, 1)
    setAll!(plan, :spectralLeakageCorrection, false)
    setAll!(plan, :sf, bSFs)
    setAll!(plan, :λ, 0)
    setAll!(plan, :tfCorrection, false)  
    c2 = reconstruct(build(plan), b)


    calibsize=hcat((calibSize.(bSFs)...)).+[6, 6, 0]*ones(Int,4)'
    fov=hcat((calibFov.(bSFs)...)).+[0.006,0.006,0]*ones(Int,4)'         
    setAll!(plan, :opParams, RecoPlan(ExplicitMultiPatchParameter))
    setAll!(plan, :tfCorrection, false)
    setAll!(plan, :gridsize, calibsize)
    setAll!(plan, :fov, fov)
    setAll!(plan, :mapping, 1:4)
    c_extr3 = reconstruct(build(plan), b)                 
    @test size(c2[1,:,:,:,1]) .+ (6,6,0) == size(c_extr3[1,:,:,:,1])               
    exportImage(joinpath(imgdir, "ExtrapolatedMultiPatch1.png"), Array(c_extr3[1,:,:,1,1]))
    @test compareImg("ExtrapolatedMultiPatch1.png")

    SM[vec([495:505;527:537]),:] .= 0.0 + 0.0im
    exportImage(joinpath(imgdir, "deadPixelSM.png"), abs.(squeeze(reshape(SM[:,52],shape(grid)...))))
    rpSM1 = getSF(bSF,freq;deadPixels=vec([495:505;527:537]))[1]
    @test rpSM1[vec([495:505;527:537]),:] != 0.0 + 0.0im
    rpSM2 = repairSM(SM,grid,Tuple.(CartesianIndices(Tuple(shape(grid)))[vec([495:505;527:537])]))
    @test rpSM1 == rpSM2
    exportImage(joinpath(imgdir, "RepairedSM.png"), abs.(squeeze(reshape(rpSM1[:,52],shape(grid)...))))


    SMtest = getSF(bSF,joinpath(datadir, "calibrations", "12extrapolated.mdf"); gridsize=[40,40,1], fov=calibFov(bSF).+[0.01,0.01,0])[1]
    bSF2 = MPIFile(joinpath(datadir, "calibrations", "12extrapolated.mdf"))
    @test calibFov(bSF2)==[0.05,0.05,0.001]
    getSF(bSF2)[1] == SMtest

end