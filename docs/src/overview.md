# Overview

We will start looking at a very basic reconstruction script
```julia
using MPIReco

fSF = MPIFile("SF_MP")
f = MPIFile("dataMP01")

c1 = reconstruction(fSF, f;
                   SNRThresh=5,
                   frames=1,
                   minFreq=80e3,
                   recChannels=1:2,
                   iterations=1,
                   spectralLeakageCorrection=true)

```
Lets go through that script step by step. First, we create handles for the system
matrix and the measurement data. Both are of the type `MPIFile` which is an abstract
type that can for instance be an `MDFFile` or a `BrukerFile`.

Using the handles to the MPI datasets we can call the `reconstruction` function
that has various variants depending on the t

## Layers
