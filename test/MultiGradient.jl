using MPIReco

# multi-gradient is a special case of multi-patch
@testset "multi-gradient in-memory reconstruction" begin
  # low gradient with patch no 3
  b = MultiMPIFile(["$datadir/MG_G1", "$datadir/MG_G2_03"])
  bSFs = MultiMPIFile(["$datadir/SF_MG_G1", "$datadir/SF_MG_G2"])
  names = names = (:color, :x, :y, :z, :time)

  c1 = reconstruction(bSFs, b;
			    SNRThresh=2, frames=1, lambd=0.003, minFreq=80e3,
			    recChannels=1:2,iterations=3, roundPatches=false)
  @test axisnames(c1) == names
  @test axisvalues(c1) == (1:1, -26.0u"mm":1.0u"mm":26.0u"mm", -26.75u"mm":1.0u"mm":26.25u"mm", 0.0u"mm":1.0u"mm":0.0u"mm", 0.0u"ms":0.6528u"ms":0.0u"ms")
  im1 = reverse(c1[1,:,:,1,1]',dims=1)
  exportImage("$imgdir/MultiGradient1.png", im1)

  #Iabs = abs.(im1)
  #Iabs[Iabs .< 0.2*maximum(Iabs)] .= 0
  #Icolored = colorview(Gray, Iabs./maximum(Iabs))
  #FileIO.save("$imgdir/MultiGradient1.png", Icolored )

  # low gradient with all 3 high gradient patches
  b = MultiMPIFile(["$datadir/MG_G1", "$datadir/MG_G2_01", "$datadir/MG_G2_02", "$datadir/MG_G2_03"])

  c2 = reconstruction(bSFs, b;
			    SNRThresh=2, frames=1, lambd=0.003, minFreq=80e3,
			    recChannels=1:2,iterations=3, roundPatches=false)
  @test axisnames(c2) == names
  # TODO wo kommt dieser Floatingpoint Fehler rein? Lässt der sich vermeiden?
  @test axisvalues(c2) == (1:1, -26.0u"mm":1.0u"mm":26.0u"mm", -26.500000000000004u"mm":1.0u"mm":30.499999999999996u"mm", 0.0u"mm":1.0u"mm":0.0u"mm", 0.0u"ms":0.6528u"ms":0.0u"ms")
  im2 = reverse(c2[1,:,:,1,1]',dims=1)
  exportImage("$imgdir/MultiGradient2.png", im2)

  #Iabs = abs.(im2)
  #Iabs[Iabs .< 0.2*maximum(Iabs)] .= 0
  #Icolored = colorview(Gray, Iabs./maximum(Iabs))
  #FileIO.save("$imgdir/MultiGradient2.png", Icolored )


  # low gradient only
  b = MultiMPIFile(["$datadir/MG_G1"])
  bSFs = MultiMPIFile(["$datadir/SF_MG_G1"])

  c3 = reconstruction(bSFs, b;
			    SNRThresh=2, frames=1, lambd=0.003, minFreq=80e3,
			    recChannels=1:2,iterations=3, roundPatches=false)
  im3 = reverse(c3[1,:,:,1,1]',dims=1)
  exportImage("$imgdir/MultiGradient3.png", im3)

  #Iabs = abs.(im3)
  #Iabs[Iabs .< 0.2*maximum(Iabs)] .= 0
  #Icolored = colorview(Gray, Iabs./maximum(Iabs))
  #FileIO.save("$imgdir/MultiGradient3.png", Icolored )
end
