RecoUtils.toDictValue(file::MPIFile) = filepath(file)
RecoUtils.toTOML(norm::AbstractRegularizationNormalization) = string(typeof(norm))