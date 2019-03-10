var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MPIReco.jl-1",
    "page": "Home",
    "title": "MPIReco.jl",
    "category": "section",
    "text": "Julia package for the reconstruction of magnetic particle imaging (MPI) data"
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "This project provides functions for the reconstruction of MPI data. The project is implemented in the programming language Julia and contains algorithms forBasic Reconstruction using a system matrix based approach\nMulti-Patch Reconstruction for data that has been acquired using a focus field sequence\nMulti-Contrast Reconstruction\nMatrix-Compression TechniquesKey features areFrequency filtering for memory efficient reconstruction. Only frequencies used during reconstructions are loaded into memory.\nDifferent solvers provided by the package RegularizedLeastSquares.jl\nHigh-level until low-level reconstruction providing maximum flexibility for the user\nSpectral leakage correction (implemented in MPIFiles.jl)"
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Start julia and open the package mode by entering ]. Then enteradd MPIRecoThis will install the packages MPIReco.jl and all its dependencies. In particular this will install the core dependencies MPIFiles and RegularizedLeastSquares."
},

{
    "location": "index.html#License-/-Terms-of-Usage-1",
    "page": "Home",
    "title": "License / Terms of Usage",
    "category": "section",
    "text": "The source code of this project is licensed under the MIT license. This implies that you are free to use, share, and adapt it. However, please give appropriate credit by citing the project."
},

{
    "location": "index.html#Contact-1",
    "page": "Home",
    "title": "Contact",
    "category": "section",
    "text": "If you have problems using the software, find mistakes, or have general questions please use the issue tracker to contact us."
},

{
    "location": "index.html#Contributors-1",
    "page": "Home",
    "title": "Contributors",
    "category": "section",
    "text": "Tobias Knopp\nMartin Möddel\nPatryk Szwargulski"
},

{
    "location": "overview.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "overview.html#Getting-Started-1",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "section",
    "text": "In order to get started we will first gather some MPI data. To this end we enter the Pkg mode in Julia (]) and execute the unit tests of MPIRecotest MPIRecoNow there will be several MPI files in the test directory. All the following examples assume that you entered the test directory and loaded MPIReco usingusing MPIReco\ncd(joinpath(dirname(pathof(MPIReco)),\"..\",\"test\"))"
},

{
    "location": "overview.html#First-Reconstruction-1",
    "page": "Getting Started",
    "title": "First Reconstruction",
    "category": "section",
    "text": "We will start looking at a very basic reconstruction scriptusing MPIReco\n\nfSF = MPIFile(\"SF_MP\")\nf = MPIFile(\"dataMP01\")\n\nc = reconstruction(fSF, f;\n                   SNRThresh=5,\n                   frames=1:10,\n                   minFreq=80e3,\n                   recChannels=1:2,\n                   iterations=1,\n                   spectralLeakageCorrection=true)\nLets go through that script step by step. First, we create handles for the system matrix and the measurement data. Both are of the type MPIFile which is an abstract type that can for instance be an MDFFile or a BrukerFile.Using the handles to the MPI datasets we can call the reconstruction function that has various variants depending on the types that are passed to it. Here, we exploit the multiple dispatch mechanism of julia. In addition to the file handles we also apply several reconstruction parameters using keyword arguments. In this case, we set the SNR threshold to 5 implying that only matrix rows with an SNR above 5 are used during reconstruction. The parameter frame decides which frame of the measured data should be reconstructed.The object c is of type ImageMeta and contains not only the reconstructed data but also several metadata such as the reconstruction parameters being used. More details on the return type are discussed in the Reconstruction Results"
},

{
    "location": "overview.html#Data-Storage-1",
    "page": "Getting Started",
    "title": "Data Storage",
    "category": "section",
    "text": "One can store the reconstruction result into an MDF file by callingsaveRecoDataMDF(\"filename.mdf\", c)In order to load the data one callsc = loaddata(\"filename.mdf\", c)We will next take a closer look at different forms of the reconstruction routine."
},

{
    "location": "basicReconstruction.html#",
    "page": "Basic Reconstruction",
    "title": "Basic Reconstruction",
    "category": "page",
    "text": ""
},

{
    "location": "basicReconstruction.html#Basic-Reconstruction-1",
    "page": "Basic Reconstruction",
    "title": "Basic Reconstruction",
    "category": "section",
    "text": "MPIReco.jl provides different reconstruction levels. All of these reconstruction routines are called reconstruction and the dispatch is done based on the input types."
},

{
    "location": "basicReconstruction.html#On-Disk-Reconstruction-1",
    "page": "Basic Reconstruction",
    "title": "On Disk Reconstruction",
    "category": "section",
    "text": "This is the highest level reconstruction. The function signature is given byfunction reconstruction(d::MDFDatasetStore, study::Study,\n                        exp::Experiment, recoParams::Dict)This reconstruction is also called an on disk reconstruction because it assumes that one has a data store (i.e. a structured folder of files) where the file location is uniquely determined by the study name and experiment number. All reconstruction parameters are passed to this method by the recoParams dictionary. On disk reconstruction has the advantage that the routine will perform reconstruction only once for a particular set of parameters. If that parameter set has already been reconstructed, the data will loaded from disk. However, the on disk reconstruction needs some experience with dataset stores to set it up correctly and is not suited for unstructured data."
},

{
    "location": "basicReconstruction.html#In-Memory-Reconstruction-1",
    "page": "Basic Reconstruction",
    "title": "In Memory Reconstruction",
    "category": "section",
    "text": "The next level is the in memory reconstruction. Its function signature readsfunction reconstruction(recoParams::Dict)This routine requires that all parameters are put into a dictionary. An overview how this dictionary looks like is given in the section Parameters.The above reconstruction method basically does two thingsPull out the location of measurement data and system matrix from the recoParams dictionary.\nPass all parameter to the low level reconstruction method in the form of keyword parameters.In turn the next level reconstruction looks like thisfunction reconstruction(bSF::Union{T,Vector{T}}, bMeas::MPIFile; kargs...)There are, however also some reconstruction methods in-between that look like thisfunction reconstruction(filenameSF::AbstractString, filenameMeas::AbstractString; kargs...)\nfunction reconstruction(filenameMeas::AbstractString; kargs...)In both cases, an MPIFile is created based on the input filename. The second version also guesses the system matrix based on what is stored within the measurement file. This usually only works, if this is executed on a system where the files are stored at exactly the same location as how they have been measured."
},

{
    "location": "basicReconstruction.html#Middle-Layer-Reconstruction-1",
    "page": "Basic Reconstruction",
    "title": "Middle Layer Reconstruction",
    "category": "section",
    "text": "The middle level reconstruction first checks, whether the dataset is a multi-patch or a single-patch file. Then it will call either reconstructionSinglePatch or reconstructionMultiPatch. Both have essentially the signaturefunction reconstructionSinglePatch(bSF::Union{T,Vector{T}}, bMeas::MPIFile;\n                                  minFreq=0, maxFreq=1.25e6, SNRThresh=-1,\n                                  maxMixingOrder=-1, numUsedFreqs=-1, sortBySNR=false, recChannels=1:numReceivers(bMeas),\n                                  bEmpty = nothing, bgFrames = 1, fgFrames = 1,\n                                  varMeanThresh = 0, minAmplification=2, kargs...) where {T<:MPIFile}"
},

{
    "location": "parameters.html#",
    "page": "Parameters",
    "title": "Parameters",
    "category": "page",
    "text": ""
},

{
    "location": "parameters.html#Parameters-1",
    "page": "Parameters",
    "title": "Parameters",
    "category": "section",
    "text": ""
},

{
    "location": "recoResults.html#",
    "page": "Results",
    "title": "Results",
    "category": "page",
    "text": ""
},

{
    "location": "recoResults.html#Reconstruction-Results-1",
    "page": "Results",
    "title": "Reconstruction Results",
    "category": "section",
    "text": "The object c is of type ImageMeta and contains not only the reconstructed data but also several metadata such as the reconstruction parameters being used. c has in total 5 dimensions. The first dimension encodes multi-spectral channels. Dimensions 2-4 encode the three spatial dimensions. The last dimension contains the number of frames being stored in c."
},

{
    "location": "multiContrast.html#",
    "page": "Multi-Contrast",
    "title": "Multi-Contrast",
    "category": "page",
    "text": ""
},

{
    "location": "multiContrast.html#Multi-Contrast-Reconstruction-1",
    "page": "Multi-Contrast",
    "title": "Multi-Contrast Reconstruction",
    "category": "section",
    "text": ""
},

{
    "location": "multiPatch.html#",
    "page": "Multi-Patch",
    "title": "Multi-Patch",
    "category": "page",
    "text": ""
},

{
    "location": "multiPatch.html#Multi-Patch-Reconstruction-1",
    "page": "Multi-Patch",
    "title": "Multi-Patch Reconstruction",
    "category": "section",
    "text": ""
},

{
    "location": "matrixCompression.html#",
    "page": "Compression",
    "title": "Compression",
    "category": "page",
    "text": ""
},

{
    "location": "matrixCompression.html#Matrix-Compression-Techniques-1",
    "page": "Compression",
    "title": "Matrix-Compression Techniques",
    "category": "section",
    "text": ""
},

]}
