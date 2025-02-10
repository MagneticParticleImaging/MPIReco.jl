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

  params = Dict{Symbol, Any}()
  params[:sf] = bSF
  params[:frames] = 1:2
  params[:alpha] = alpha
  params[:choosePeak] = choosePeak
  params[:windowType] = windowType
  params[:higherHarmonic] = choosePeak
  params[:lambda] = lambda
  params[:iterations] = iterations # reconstruction parameter
  params[:SNRThresh] = 10
  params[:minFreq] = 80e3
  params[:recChannels] = 1:3
  params[:λ] = 0.0f0

  c = reconstruct("MultiPatchMotion", bMeas; params...)
  exportImage(joinpath(imgdir, "Motion1.png"), maximum(Array(c[1,:,:,:,1]),dims=3)[:,:,1])
  @test compareImg("Motion1.png")

  c = reconstruct("MultiPatchMotion", bMeas; params..., bgParams = SimpleExternalBackgroundCorrectionParameters(;emptyMeas = bBG, bgFrames = bgFrames))
  exportImage(joinpath(imgdir, "Motion2.png"), maximum(Array(c[1,:,:,:,1]),dims=3)[:,:,1])
  @test compareImg("Motion2.png")

end
