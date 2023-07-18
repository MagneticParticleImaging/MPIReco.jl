using MPIReco

@testset "multi-patch in-memory reconstruction" begin
  bSF = MultiMPIFile([joinpath(datadir, "calibrations", "12.mdf")])
  dirs = ["1.mdf", "2.mdf", "3.mdf", "4.mdf"]
  b = MultiMPIFile(joinpath.(datadir, "measurements", "20211226_203916_MultiPatch", dirs))
  names = (:color, :x, :y, :z, :time)
  values1 = (1:1,
	    -27.5u"mm":1.25u"mm":27.5u"mm",
	    -27.5u"mm":1.25u"mm":27.5u"mm",
	    0.0u"mm":1.0u"mm":0.0u"mm",
	    0.0u"ms":0.6528u"ms":0.0u"ms")
  values2 = (values1[1], -19.375u"mm":1.25u"mm":19.375u"mm",
	     -19.375u"mm":1.25u"mm":19.375u"mm", values1[4:5]...)

  # basic multi-patch reconstruction
  c1 = reconstruction(bSF, b;
			    SNRThresh=5, frames=1, minFreq=80e3,
			    recChannels=1:2,iterations=1,
			    spectralLeakageCorrection=false)
  @test axisnames(c1) == names
  @test axisvalues(c1) == values1
  exportImage(joinpath(imgdir, "MultiPatch1.png"), arraydata(c1[1,:,:,1,1]))
  @test compareImg("MultiPatch1.png")

  # TODO test description
  c2 = reconstruction(bSF, b;
			    SNRThresh=5, frames=1, minFreq=80e3,
			    recChannels=1:2,iterations=1, roundPatches=true,
			    spectralLeakageCorrection=false)
  @test axisnames(c2) == names
  @test axisvalues(c2) == values1
  exportImage(joinpath(imgdir, "MultiPatch2.png"), arraydata(c2[1,:,:,1,1]))
  @test compareImg("MultiPatch2.png")

  # multi-patch reconstruction using multiple system matrices
  dirs = ["8.mdf", "9.mdf", "10.mdf", "11.mdf"]
  bSFs = MultiMPIFile(joinpath.(datadir, "calibrations", dirs))
  c3 = reconstruction(bSFs, b;
			    SNRThresh=5, frames=1, minFreq=80e3,
			    recChannels=1:2,iterations=1,
			    spectralLeakageCorrection=false)
  @test axisnames(c3) == names
  @test axisvalues(c3) == values2
  exportImage(joinpath(imgdir, "MultiPatch3.png"), arraydata(c3[1,:,:,1,1]))
  @test compareImg("MultiPatch3.png")

  # flexible multi-patch reconstruction
  dirs = ["8.mdf", "9.mdf", "10.mdf", "11.mdf"]
  bSFs = MultiMPIFile(joinpath.(datadir, "calibrations", dirs))
  mapping = [1,2,3,4]
  freq = filterFrequencies(bSFs, SNRThresh=5, minFreq=80e3)
  S = [getSF(SF,freq,nothing,"kaczmarz", bgcorrection=false)[1] for SF in bSFs]
  SFGridCenter = zeros(3,4)
  FFPos = zeros(3,4)
  FFPos[:,1] = [-0.008, 0.008, 0.0]
  FFPos[:,2] = [-0.008, -0.008, 0.0]
  FFPos[:,3] = [0.008, 0.008, 0.0]
  FFPos[:,4] = [0.008, -0.008, 0.0]
  c4 = reconstruction(bSFs, b;
			    SNRThresh=5, frames=1, minFreq=80e3,
			    recChannels=1:2,iterations=1,
			    spectralLeakageCorrection=false, mapping=mapping,
			    systemMatrices = S, SFGridCenter=SFGridCenter,
			    FFPos=FFPos, FFPosSF=FFPos)
  @test axisnames(c4) == names
  @test axisvalues(c4) == values2
  exportImage(joinpath(imgdir, "MultiPatch4.png"), arraydata(c4[1,:,:,1,1]))
  @test compareImg("MultiPatch4.png")
  # TODO the last test shows odd results
end
