# Background Estimation Example

This folder contains example code for applying a joint reconstruction of
background and particle concentration of MPI data. The method uses a background
dictionary that is setup using the background scans measured during system
matrix acquisition.

## Installation

In order to use this code one first has to download Julia (version 1.x) and install
`MPIReco` by executing

```julia
using Pkg
Pkg.add("MPIReco")
```

Load the package by entering
```julia
using MPIReco
```

You can then switch to the directory of this example by entering
```julia
dir = joinpath(dirname(pathof(MPIReco)), "..","examples","BGEstimation")
cd(dir)
```

## Execution
After installation and switching to the example directory, the example code can be
executed by entering

```julia
include("example.jl")
```

This will first download all data (about 134 MB) and then perform a reconstruction.
Parameters of the reconstruction are documented in the Julia script and can be
changed. After the reconstruction is done, the script will open a PyPlot window
and show the reconstruction result. 
