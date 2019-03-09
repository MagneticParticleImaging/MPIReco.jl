var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#MPI-Reco-1",
    "page": "Home",
    "title": "MPI Reco",
    "category": "section",
    "text": "Julia package for reconstruction of magnetic particle imaging (MPI) data"
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "This project provides functions for the reconstruction of MPI data. The project is implemented in the programming language Julia and contains algorithms forregular 1D / 2D / 3D image reconstruction using a system matrix based approach\nmulti-patch and multi-gradient reconstruction for data the has been acquired using a focus field sequence\nmulti-colored image reconstruction\nmatrix-compression techniquesKey features arefrequency filtering for memory efficient reconstruction. Only frequencies used during reconstructions are loaded into memory.\nspectral leakage correction (implemented in MPIFiles.jl)"
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
    "text": "Florian Griese\nNadine Gdaniec\nTobias Knopp\nMartin Möddel\nPatryk Szwargulski"
},

{
    "location": "usage.html#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "usage.html#Installation-and-Usage-1",
    "page": "Installation",
    "title": "Installation and Usage",
    "category": "section",
    "text": ""
},

{
    "location": "usage.html#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "The MPIReco project is a package written in the programming language Julia. In order to install MPIReco you first have to install Julia in version 0.7 (version 1.0 not yet supported). Then open Julia and enter ] to open the package mode. Then enterdev https://github.com/MagneticParticleImaging/MPIFiles.jl.git\ndev https://github.com/tknopp/RegularizedLeastSquares.jl.git\ndev https://github.com/MagneticParticleImaging/MPIReco.jl.gitwhich will install the package including its dependencies MPIFiles and RegularizedLeastSquares. Then enterusing MPIRecoto load the packageIn order to obtain the example datasets you have to execute the unit tests which can be done by enteringtest MPIRecowithin the package mode of Julia."
},

{
    "location": "usage.html#Getting-Started-1",
    "page": "Installation",
    "title": "Getting Started",
    "category": "section",
    "text": "In order to get started you can executeusing Pkg\ninclude(joinpath(dirname(pathof(MPIReco)),\"..\",\"test\",\"MultiPatch.jl\"))which will perform exemplary reconstructions and plot the results."
},

{
    "location": "overview.html#",
    "page": "Overview",
    "title": "Overview",
    "category": "page",
    "text": ""
},

{
    "location": "overview.html#Overview-1",
    "page": "Overview",
    "title": "Overview",
    "category": "section",
    "text": "We will start looking at a very basic reconstruction scriptusing MPIReco\n\nfSF = MPIFile(\"SF_MP\")\nf = MPIFile(\"dataMP01\")\n\nc = reconstruction(fSF, f;\n                   SNRThresh=5,\n                   frames=1:10,\n                   minFreq=80e3,\n                   recChannels=1:2,\n                   iterations=1,\n                   spectralLeakageCorrection=true)\nLets go through that script step by step. First, we create handles for the system matrix and the measurement data. Both are of the type MPIFile which is an abstract type that can for instance be an MDFFile or a BrukerFile.Using the handles to the MPI datasets we can call the reconstruction function that has various variants depending on the types that are passed to it. Here, we exploit the multiple dispatch mechanism of julia. In addition to the file handles we also apply several reconstruction parameters using keyword arguments. In this case, we set the SNR threshold to 5 implying that only matrix rows with an SNR above 5 are used during reconstruction. The parameter frame decides which frame of the measured data should be reconstructed.The object c is of type ImageMeta and contains not only the reconstructed data but also several metadata such as the reconstruction parameters being used. c has in total 5 dimensions. The first dimension encodes multi-spectral channels. Dimensions 2-4 encode the three spatial dimensions. The last dimension contains the number of frames being stored in c."
},

{
    "location": "overview.html#Data-Storage-1",
    "page": "Overview",
    "title": "Data Storage",
    "category": "section",
    "text": "One can store the reconstruction result into an MDF file by callingsaveRecoDataMDF(\"filename.mdf\", c)In order to load the data one callsc = loaddata(\"filename.mdf\", c)"
},

{
    "location": "multiPatch.html#",
    "page": "Multi-patch",
    "title": "Multi-patch",
    "category": "page",
    "text": ""
},

]}
