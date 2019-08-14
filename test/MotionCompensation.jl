using MPIReco

@testset "motion compensation reconstruction" begin

  ################
  ## Parameters ##
  ################
  alpha = 0.4 

  # Choose the peak number according to the expected approximate value. Be careful, if there are peaks not corresponding to object motion, division by choosePeak can be false.
  # For experiments with different periodic motions it is used to pick one motion frequency (then don't divide by choosePeak)

  choosePeak = 4  

  windowType = 1  # 1: Hann, 2: FT1A05, 3: Rectangle

  # Measurement data
  datadirMeas = "./motionComp/"
  b = MPIFile(datadirMeas*"measFast.mdf") # high frequency
  bBG = MPIFile(datadirMeas*"measBG.mdf") # background measurement

  # System matrices
  datadirSF = "./motionComp/"
  SFall = ["SF1Small.mdf","SF2Small.mdf","SF3Small.mdf","SF4Small.mdf"]
  bSF = MultiMPIFile(datadirSF.*SFall)

  # Frequency selection # TODO: reconstructionPeriodicMotion basteln, die das tut
  freq = filterFrequencies(bSF,SNRThresh=2,minFreq=80e3,recChannels=[1,2,3])

  # Background frames
  frBG = [1:200,201:400, 401:600 ,601:800]

  # Selected frames
  recoFrame = 1

  ####################
  ## Reconstruction ##
  ####################

  # Reconstruction parameters
  lambda = 0.01 
  iterations = 1

  c = reconstructionPeriodicMotion(b, bSF, bBG,
		frBG, choosePeak, alpha, freq, recoFrame,
		lambda=lambda, iterations=iterations,
		windowType=windowType, higherHarmonic=choosePeak)

  exportImage("./img/MotionComp.png", maximum(arraydata(c[:,:,:,1]),dims=3)[:,:,1])

end
