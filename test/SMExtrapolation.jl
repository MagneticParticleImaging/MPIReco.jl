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
    c1 = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
    recChannels=1:2, iterations=1, spectralLeakageCorrection=true, SMextrapolation=nothing)

    c_extr = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
    recChannels=1:2, iterations=1, spectralLeakageCorrection=true, gridsize=[36,36,1], fov=calibFov(bSF).+[0.006,0.006,0]) 
    @test size(c1[1,:,:,:,1]) .+ (4,4,0) == size(c_extr[1,:,:,:,1])
    exportImage(joinpath(imgdir, "Extrapolated1.png"), arraydata(c_extr[1,:,:,1,1]))

    c_extr2 = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
    recChannels=1:2, iterations=1, spectralLeakageCorrection=true, fov=calibFov(bSF).+[0.006,0.006,0])
    @test size(c1[1,:,:,:,1]) == size(c_extr2[1,:,:,:,1])
    exportImage(joinpath(imgdir, "Extrapolated2.png"), arraydata(c_extr2[1,:,:,1,1]))

    SFdirs = ["8.mdf", "9.mdf", "10.mdf", "11.mdf"]
    bSFs = MultiMPIFile(joinpath.(datadir, "calibrations", SFdirs))
    bdirs = ["1.mdf", "2.mdf", "3.mdf", "4.mdf"]
    b = MultiMPIFile(joinpath.(datadir, "measurements", "20211226_203916_MultiPatch", bdirs))
    c2 = reconstruction(bSFs, b;
                  SNRThresh=5, frames=1, minFreq=80e3,
                  recChannels=1:2,iterations=1,
                  spectralLeakageCorrection=false
                  )
    calibsize=hcat((calibSize.(bSFs)...)).+[6, 6, 0]*ones(Int,4)'
    fov=hcat((calibFov.(bSFs)...)).+[0.006,0.006,0]*ones(Int,4)'         
    c_extr3 = reconstruction(bSFs, b;
                  SNRThresh=5, frames=1, minFreq=80e3,
                  recChannels=1:2,iterations=1,
                  spectralLeakageCorrection=false,
                  gridsize=calibsize,fov=fov)                  
    @test size(c2[1,:,:,:,1]) .+ (6,6,0) == size(c_extr3[1,:,:,:,1])               
    exportImage(joinpath(imgdir, "ExtrapolatedMultiPatch1.png"), arraydata(c_extr3[1,:,:,1,1]))

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