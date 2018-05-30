# Overview

We will start looking at a very basic reconstruction script
```julia
using MPIReco

fSF = MPIFile("SF_MP")
f = MPIFile("dataMP01")

c = reconstruction(fSF, f;
                   SNRThresh=5,
                   frames=1:10,
                   minFreq=80e3,
                   recChannels=1:2,
                   iterations=1,
                   spectralLeakageCorrection=true)

```
Lets go through that script step by step. First, we create handles for the system
matrix and the measurement data. Both are of the type `MPIFile` which is an abstract
type that can for instance be an `MDFFile` or a `BrukerFile`.

Using the handles to the MPI datasets we can call the `reconstruction` function
that has various variants depending on the types that are passed to it. Here, we
exploit the multiple dispatch mechanism of julia. In addition to the file handles
we also apply several reconstruction parameters using keyword arguments. In this case,
we set the SNR threshold to 5 implying that only matrix rows with an SNR above 5 are used
during reconstruction. The parameter frame decides which frame of the measured data
should be reconstructed.

The object c is of type `ImageMeta` and contains not only the reconstructed data
but also several metadata such as the reconstruction parameters being used.
c has in total 5 dimensions. The first dimension encodes multi-spectral channels.
Dimensions 2-4 encode the three spatial dimensions. The last dimension contains the
number of frames being stored in c.

## Data Storage

One can store the reconstruction result into an MDF file by calling
```julia
saveRecoDataMDF("filename.mdf", c)
```
In order to load the data one calls
```julia
c = loaddata("filename.mdf", c)
```
