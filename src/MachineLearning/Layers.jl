### ConvAsymPadding Layer ###

struct ConvAsymPadding{T, N}
    c::T
    pad::NTuple{N,Int}
end
function ConvAsymPadding(k::NTuple{N,Integer}, ch::Pair{<:Integer,<:Integer}, σ = identity;
    init = Flux.glorot_uniform,  stride = 1, pad = 0, dilation = 1, bias=true) where N
    length(pad) < 2N || all(i -> length(unique(pad[i:i+1])) == 1, 1:2:2N) == 1&& return Conv(k, ch , σ, init=init, stride=stride, pad=pad, dilation=dilation, bias=bias)

    pad_manual = Tuple(map(i -> abs(pad[i] - pad[i+1]), 1:2:2N))
    pad_auto = Tuple(map(i -> minimum(pad[i:i+1]), 1:2:2N))
    return ConvAsymPadding(Conv(k, ch, σ, init=init, stride=stride, pad=pad_auto, dilation=dilation), pad_manual)
end
function (c::ConvAsymPadding)(x::AbstractArray)
    # Maybe there are faster ways to do this as well...
    padding = similar(x, c.pad..., size(x)[end-1:end]...)
    fill!(padding, 0)
    c.c(cat(x, padding, dims=1:length(c.pad)))
end

Flux.@functor ConvAsymPadding


### Cropping Layer ###

struct Crop{N}
    crop::NTuple{N,Int}
end

function (c::Crop{4})(x::AbstractArray)
  return x[(1+c.crop[1]):(end-c.crop[2]), (1+c.crop[3]):(end-c.crop[4]),:,:]
end

function (c::Crop{6})(x::AbstractArray)
  return x[(1+c.crop[1]):(end-c.crop[2]), (1+c.crop[3]):(end-c.crop[4]), (1+c.crop[5]):(end-c.crop[6]),:,:]
end

Flux.@functor Crop

### BatchNormWrap ###

expand_dims(x,n::Int) = reshape(x,ones(Int64,n)...,size(x)...)
	
function _random_normal(shape...)
  return Float32.(rand(Normal(0.f0,0.02f0),shape...)) |> device
end
	
function BatchNormWrap(out_ch)
    Chain(x->expand_dims(x,3),
	  BatchNorm(out_ch),
	  x->reshape(x, size(x)[4:end]))
end

### UNetConvBlock ###

UNetConvBlock(in_chs, out_chs; kernel = (3,3,3), pad = (1,1,1), stride=(1,1,1)) =
  Chain(Conv(kernel, in_chs=>out_chs, pad = pad, stride=stride; init=_random_normal),
	BatchNormWrap(out_chs),
	x->leakyrelu.(x,0.02f0))