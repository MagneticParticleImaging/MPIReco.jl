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
  datadirMeas = joinpath(datadir, "mdf")
  bMeas = MPIFile(joinpath(datadirMeas, "measFast.mdf")) # high frequency
  bBG = MPIFile(joinpath(datadirMeas, "measBG.mdf")) # background measurement

  # System matrices
  datadirSF = joinpath(datadir, "mdf")
  SFall = ["SF1Small.mdf","SF2Small.mdf","SF3Small.mdf","SF4Small.mdf"]
  bSF = MultiMPIFile(joinpath.(datadirSF, SFall))

  # Background frames
  bgFrames = [1:200, 201:400, 401:600 ,601:800]

  # Selected frames
  recoFrame = 1

  ####################
  ## Reconstruction ##
  ####################

  # Reconstruction parameters
  lambda = 0.01
  iterations = 1

  c = reconstruction(bSF, bMeas, periodicMotionCorrection = true,
                                bEmpty = bBG, bgFrames = bgFrames, #background measurement
                                alpha = alpha, choosePeak = choosePeak, recoFrame = recoFrame,
                                windowType=windowType, higherHarmonic=choosePeak,
                                lambda=lambda, iterations=iterations, # reconstruction parameter
                                SNRThresh=10, minFreq=80e3, recChannels=[1,2,3]) #frequency selection

  exportImage(joinpath(imgdir, "MotionComp.png"), maximum(arraydata(c[1,:,:,:,1]),dims=3)[:,:,1])

end
