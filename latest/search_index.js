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
    "text": "The source code of this project is licensed under the MIT license. This implies that you are free to use, share, and adapt it. However, you must give appropriate credit by citing the project."
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
    "text": "The MPIReco project is a package written in the programming language Julia. In order to install MPIReco you first have to install Julia in version 0.6. Then open Julia and enterPkg.clone(\"https://github.com/MagneticParticleImaging/MPIReco.jl.git\")which will install the package. Then enterusing MPIRecowhich will install the dependencies MPIFiles and LinearSolver. Further dependencies have already been installed during the clone of the package.In order to obtain the example datasets you have to execute the unit tests which can be done by enteringPkg.test(\"MPIReco\")"
},

{
    "location": "usage.html#Getting-Started-1",
    "page": "Installation",
    "title": "Getting Started",
    "category": "section",
    "text": "In order to get started you can execute oninclude(Pkg.dir(\"MPIReco\",\"test/MultiPatch.jl\"))which will perform exemplary reconstructions and plot the results."
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
    "text": "We will start looking at a very basic reconstruction scriptusing MPIReco\n\nfSF = MPIFile(\"SF_MP\")\nf = MPIFile(\"dataMP01\")\n\nc1 = reconstruction(fSF, f;\n                   SNRThresh=5,\n                   frames=1,\n                   minFreq=80e3,\n                   recChannels=1:2,\n                   iterations=1,\n                   spectralLeakageCorrection=true)\nLets go through that script step by step. First, we create handles for the system matrix and the measurement data. Both are of the type MPIFile which is an abstract type that can for instance be an MDFFile or a BrukerFile.Using the handles to the MPI datasets we can call the reconstruction function that has various variants depending on the t"
},

{
    "location": "overview.html#Layers-1",
    "page": "Overview",
    "title": "Layers",
    "category": "section",
    "text": ""
},

{
    "location": "multiPatch.html#",
    "page": "Multi-patch",
    "title": "Multi-patch",
    "category": "page",
    "text": ""
},

]}
