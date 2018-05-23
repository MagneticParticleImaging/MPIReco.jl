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
    "text": "Initiative for open magnetic particle imaging data"
},

{
    "location": "index.html#Introduction-1",
    "page": "Home",
    "title": "Introduction",
    "category": "section",
    "text": "Magnetic particle imaging (MPI) is a tomographic imaging technique that allows to determine the spatial distribution of magnetic nanoparticles (MNPs). The aim of this project is to provide MPI data that are freely available for the research community. To facilitate reproducible research all available information including a detailed description of the used phantoms and the applied measurement sequences are provided.All MPI data is stored in the MPI Data Format (MDF). The MDF provides a common data format for the storage of MPI raw data, calibration data, and reconstruction data.This projects comes with a first set of datasets covering the most interesting cases to get started with MPI. New datasets can be easily added and researcher worldwide are invited to contribute own datasets."
},

{
    "location": "index.html#License-/-Terms-of-Usage-1",
    "page": "Home",
    "title": "License / Terms of Usage",
    "category": "section",
    "text": "The source code of this project is licensed under the MIT license. This implies that you are free to use, share, and adapt the data. However, you must give appropriate credit by citing the project."
},

{
    "location": "index.html#Contact-1",
    "page": "Home",
    "title": "Contact",
    "category": "section",
    "text": "If you have problems using the datasets, find mistakes, or have general questions please use the issue tracker to contact us. You can also contact Tobias Knopp directly."
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
    "page": "Usage",
    "title": "Usage",
    "category": "page",
    "text": ""
},

{
    "location": "usage.html#Installation-and-Usage-1",
    "page": "Usage",
    "title": "Installation and Usage",
    "category": "section",
    "text": ""
},

{
    "location": "usage.html#Installation-1",
    "page": "Usage",
    "title": "Installation",
    "category": "section",
    "text": "The MPIReco project is a package written in the programming language Julia. In order to install MPIReco you first have to install Julia in version 0.6. Then open Julia and enterPkg.clone(\"https://github.com/MagneticParticleImaging/MPIReco.jl.git\")which will install the package. Then enterusing MPIRecowhich will install the dependencies MPIFiles and LinearSolver. Further dependencies have already been installed during the clone of the package.In order to obtain the example datasets you have to execute the unit tests which can be done by enteringPkg.test(\"MPIReco\")"
},

{
    "location": "usage.html#Getting-Started-1",
    "page": "Usage",
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

]}
