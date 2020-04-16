# System Matrix Recovery Example

This folder contains example code for the Compressed Sensing based recovery of System Matrices
from undersampled measurements.

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

You can the switch to the directory of this example by entering
```julia
dir = joinpath(dirname(pathof(MPIReco)), "..","examples","SMRecovery")
cd(dir)
```

## Execution
After installation and switching to the example directory, the example code can be
executed by entering

```julia
include("example.jl")
```

This will first download all data (about 31 GB). Then a subset of frequency components is selected
and matrix recovery is performed using three different techniques (CS, CSLR and CSFR). 
Parameters of the reconstruction are documented in the Julia script and can be
changed. After the reconstruction is done, the script will open a PyPlot window
and show slices through the recovered frequency components.
