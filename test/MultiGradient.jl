using MPIReco
using Test
using Winston

# multi-gradient is a special case of multi-patch
@testset "multi-gradient in-memory reconstruction" begin
  # low gradient with patch no 3
  b = MultiMPIFile(["dataMG_G1", "dataMG_G2_03"])
  bSFs = MultiMPIFile(["SF_MG_G1", "SF_MG_G2"])
  names = names = (:color, :x, :y, :z, :time)
  Winston.colormap("Grays")

  c1 = reconstruction(bSFs, b;
			    SNRThresh=2, frames=1, lambd=0.003, minFreq=80e3,
			    recChannels=1:2,iterations=3, roundPatches=false)
  @test axisnames(c1) == names
  @test axisvalues(c1) == (1:1, -26.0u"mm":1.0u"mm":26.0u"mm", -26.75u"mm":1.0u"mm":26.25u"mm", 0.0u"mm":1.0u"mm":0.0u"mm", 0.0u"ms":0.6528u"ms":0.0u"ms")
  im1 = reverse(c1[1,:,:,1,1]',dims=1)
  p = imagesc(data(data(im1)), (minimum(im1),maximum(im1)))
  #display(p)
  
  # low gradient with all 3 high gradient patches
  b = MultiMPIFile(["dataMG_G1", "dataMG_G2_01", "dataMG_G2_02", "dataMG_G2_03"])
  
  c2 = reconstruction(bSFs, b;
			    SNRThresh=2, frames=1, lambd=0.003, minFreq=80e3,
			    recChannels=1:2,iterations=3, roundPatches=false)
  @test axisnames(c2) == names
  # TODO wo kommt dieser Floatingpoint Fehler rein? Lässt der sich vermeiden?
  @test axisvalues(c2) == (1:1, -26.0u"mm":1.0u"mm":26.0u"mm", -26.500000000000004u"mm":1.0u"mm":30.499999999999996u"mm", 0.0u"mm":1.0u"mm":0.0u"mm", 0.0u"ms":0.6528u"ms":0.0u"ms")
  im2 = reverse(c2[1,:,:,1,1]',dims=1)
  p = imagesc(data(data(im2)), (minimum(im2),maximum(im2)))
  #display(p)
end
