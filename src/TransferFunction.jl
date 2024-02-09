
export estimateTransferFunction, correctTransferFunction

# Estimate the transfer function of the receive chain using linear regression
function estimateTransferFunction(SMeas, SModel)
  a = zeros(ComplexF64, size(SMeas,3), size(SMeas,4))
  for k in CartesianIndices(a)
    a[k] = sum(conj(vec(SModel[:,:,k])) .* vec(SMeas[:,:,k])) / norm(vec(SModel[:,:,k])).^2
  end
  return a
end

function correctTransferFunction(SMeas, SModel)
  a = estimateTransferFunction(SMeas, SModel)
  SModel = reshape(a,1,1,size(SModel,3),size(SModel,4)) .* SModel
  SModel[:,:,1,:] .= 1
  return SModel
end


function estimateTransferFunction(SMeas, SModel, tfMeasured, shift, phi)

  transfer_signs = [1im,1im,1im]
#  transfer_signs = [1,1,1]

  # shift correction
  N_ = (size(SModel,3)-1)*2
  shiftTerm = reshape( exp.(2*pi*im*phi .-2*pi*im.*(0:(size(SModel,3)-1))/N_*shift ),:,1)

  # here we to the scaling using mixing factor mx=my=3
  fac = maximum(abs.(SMeas[:, :, (3-1)*16+(3-1)*17+1, :]),dims=(1,2)) ./ 
      maximum(abs.(SModel[:, :, (3-1)*16+(3-1)*17+1, :]),dims=(1,2))

  # build the TF
  a_ = fac[1,:,:] .* (((tfMeasured .* reshape(transfer_signs,1,:))  ) .* shiftTerm) 

  return a_
end

function correctTransferFunction(SMeas, SModel, tfMeasured, shift, phi)
  a = estimateTransferFunction(SMeas, SModel,tfMeasured, shift, phi)
  SModel = reshape(a,1,1,size(SModel,3),size(SModel,4)) .* SModel
  SModel[:,:,1,:] .= 1
  return SModel
end
