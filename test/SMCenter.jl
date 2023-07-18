using MPIReco

@testset "system matrix center estimation" begin
  #dirs = ["8.mdf", "9.mdf", "10.mdf", "11.mdf"]   8.mdf is wrong. BG is subtracted twice
  dirs = ["9.mdf", "10.mdf", "11.mdf"]
  SFs = joinpath.(datadir, "calibrations", dirs)
  #center = [[9.0,23.0,1.0],[9.0,10.0,1.0],[22.0,23.0,1.0],[22.0,10.0,1.0]]
  center = [[9.0,10.0,1.0],[22.0,23.0,1.0],[22.0,10.0,1.0]]

  for (l,SF) in enumerate(SFs)
    bSF = MPIFile(SF)
    a = findCenterOfDfFov(bSF)
    @test a == center[l]
    S = getSF(bSF,2,0,0,2)
    exportImage(joinpath(imgdir, "Center1.png"), abs.(S)[:,:,1])
    #@test compareImg("Center1.png")
    #p = imagesc(abs.(S)[:,:,1])
    #add(p, Points([a[2]-1],[a[1]-1],color="yellow",lw=5))
    #savefig(p, "$imgdir/Center1.png")
  end

  for (l,SF) in enumerate(SFs)
    bSF = MPIFile(SF)
    a = findCenterOfDfFov(bSF)
    @test a == center[l]
    S = getSF(bSF,5,6,0,1)
    exportImage(joinpath(imgdir, "Center2.png"), abs.(S)[:,:,1])
    #@test compareImg("Center2.png")
    #p = imagesc(abs.(S)[:,:,1])
    #add(p, Points([a[2]-1],[a[1]-1],color="yellow",lw=5))
    #savefig(p, "$imgdir/Center2.png")
  end
end
