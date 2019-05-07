# Matrix-Compression Techniques

The reconstruction can be accelerated by applying matrix compression. To this end, the system matrix `S` is transformed into a different domain by applying a  basis transformation on the rows of the system matrix. In `MPIReco.jl`, matrix compression can be enabled by specifying `sparseTrafo` which can be `"DCT-IV"` or `"FFT"`. 

The transformations can be restricted to the drive-field field-of-view by setting `useDFFoV = true`. The compression factor that controls how many coefficients are dropped after application of the transformation is controlled by the parameter `redFactor`. For instance a reduction factor of `redFactor = 0.01` will drop 99 % of the data.
