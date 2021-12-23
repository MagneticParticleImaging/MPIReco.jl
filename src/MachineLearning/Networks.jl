
### make U-Net ###

function make_model_unet(N)
    layers = Any[]
    maxDepth = minimum(round.(Int,log2.(N)).-1)
    numLayers = min(maxDepth, 3)
    
    H = zeros(Int, numLayers, length(N))
    H[1,:] .= N
    for l=1:numLayers-1
      H[l+1,:] .= ceil.(Int, vec(H[l,:])./2)
    end
    needCrop = Int.(isodd.(H))
    
    interp = :nearest #trilinear

    inChan = 1
    outChan = 64
    for l=1:numLayers
      push!(layers,UNetConvBlock(inChan,outChan,kernel=(3,3,3), stride=(2,2,2), pad=1))
      push!(layers,UNetConvBlock(outChan,outChan,kernel=(3,3,3), pad=1))
      inChan = outChan
      outChan *= 2
    end
    outChan = inChan
    push!(layers, Upsample(interp,scale=(2,2,2)))
    if sum(needCrop[:,end])>0
      push!(layers, Crop( (0,needCrop[end,1],0,needCrop[end,2],0,needCrop[end,3]) ) )
    end
    push!(layers, UNetConvBlock(inChan,outChan, kernel=(3,3,3), pad=1))
    push!(layers, UNetConvBlock(inChan, outChan, kernel=(1,1,1), pad=0))
    outChan = outChan ÷ 2
    for l=1:(numLayers-1)
      push!(layers, Upsample(interp,scale=(2,2,2))) 
      if sum(needCrop[:,end-l])>0
        push!(layers, Crop( (0,needCrop[end-l,1],0,needCrop[end-l,2],0,needCrop[end-l,3]) ) )
      end  
      push!(layers, UNetConvBlock(inChan,outChan,kernel=(3,3,3), pad=1))
      push!(layers, UNetConvBlock(outChan,outChan,kernel=(1,1,1), pad=0))
      inChan = outChan
      outChan = outChan ÷ 2
    end
    push!(layers, Conv((1, 1, 1), 64=>1, pad = 0, relu ;init=_random_normal))
  
    return Chain(layers...)
  end