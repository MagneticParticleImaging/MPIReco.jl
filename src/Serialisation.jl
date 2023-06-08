function RecoUtils.toDict!(dict, reg::AbstractRegularization) 
  dict["type"] = string(typeof(reg))
  for field in fieldnames(typeof(reg))
    dict[string(field)] = toDictValue(getproperty(reg, field))
  end
  return dict
end

RecoUtils.toDictValue(file::MPIFile) = filepath(file)
RecoUtils.toDictValue(norm::AbstractRegularizationNormalization) = string(typeof(norm))