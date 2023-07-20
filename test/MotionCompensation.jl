using MPIReco

@testset "motion compensation reconstruction" begin

  ################
  ## Parameters ##
  ################
  alpha = 1.2

  # Choose the peak number according to the expected approximate value. Be careful, if there are peaks not corresponding to object motion, division by choosePeak can be false.
  # For experiments with different periodic motions it is used to pick one motion frequency (then don't divide by choosePeak)

  choosePeak = 4

  windowType = 1  # 1: Hann, 2: FT1A05, 3: Rectangle

  # Measurement data
  datadirMeas = joinpath(datadir, "measurements", "20211226_204336_Rotation")
  bMeas = MPIFile(joinpath(datadirMeas, "1.mdf")) # high frequency
  bBG = MPIFile(joinpath(datadirMeas, "2.mdf")) # background measurement

  # System matrices
  datadirSF = joinpath(datadir, "calibrations")
  SFall = ["3.mdf","4.mdf","5.mdf","6.mdf"]
  bSF = MultiMPIFile(joinpath.(datadirSF, SFall))

  # Background frames
  #bgFrames = vcat(collect.([1:200, 201:400, 401:600 ,601:800])...)
  bgFrames = 1:1 # TODO improve bg correction for multi patch, bgFrames used to be ignored

  ####################
  ## Reconstruction ##
  ####################

  # Reconstruction parameters
  lambda = 0.01
  iterations = 1

  plan = getPlan("Motion")
  setAll!(plan, :sf, bSF)
  setAll!(plan, :frames, 1:1)
  setAll!(plan, :alpha, alpha)
  setAll!(plan, :choosePeak, choosePeak)
  setAll!(plan, :windowType, windowType)
  setAll!(plan, :higherHarmonic, choosePeak)
  setAll!(plan, :lambda, lambda)
  setAll!(plan, :iterations, iterations) # reconstruction parameter
  setAll!(plan, :SNRThresh, 10)
  setAll!(plan, :minFreq, 80e3)
  setAll!(plan, :recChannels, 1:3)
  setAll!(plan, :λ, 0.0f0)

  c = reconstruct(build(plan), bMeas)
  exportImage(joinpath(imgdir, "Motion1.png"), maximum(arraydata(c[1,:,:,:,1]),dims=3)[:,:,1])
  @test compareImg("Motion1.png")

  setAll!(plan, :bgParams, SimpleExternalBackgroundCorrectionParameters(;emptyMeas = bBG, bgFrames = bgFrames))
  c = reconstruct(build(plan), bMeas)
  exportImage(joinpath(imgdir, "Motion2.png"), maximum(arraydata(c[1,:,:,:,1]),dims=3)[:,:,1])
  @test compareImg("Motion2.png")


end
