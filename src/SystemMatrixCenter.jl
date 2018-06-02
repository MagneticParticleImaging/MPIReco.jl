export findCenterOfDfFov

"""This function reads the component belonging to the given mixing-factors und chanal from the given MPIfile.\n
In the case of a 2D systemfunction, the third dimension is added."""
function readSF(bSF,mx::Int,my::Int,mz::Int,recChan::Int)
  maxFreq  = rxNumFrequencies(bSF)
	freq = mixFactorToFreq(bSF, mx, my, mz)
	freq = clamp(freq,0,maxFreq-1)
	k = freq + (recChan-1)*maxFreq
  println("k = $k")
	SF, grid = getSF(bSF,[k+1],returnasmatrix = true, bgcorrection=true)
  N = calibSize(bSF)
  SF = reshape(SF,N[1],N[2],N[3])

  return SF
end

"""This function calculates the center of the DfFov of the given MPIFile bSF"""
function findCenterOfDfFov(bSF::MPIFile)
	S1 = readSF(bSF,0,2,0,1) #xdir
  S2 = readSF(bSF,2,0,0,2) #ydir
  S3 = readSF(bSF,0,0,2,3) #zdir

  center = ones(3)
  if size(S1,1) > 1
    u = floor.(Int, centerOfMass(abs.(S1)))
    println(u)
    y_0 = real(S1[u[1]-1,u[2],u[3]])
    y_1 = real(S1[u[1]+1,u[2],u[3]])
    xdiff = -y_0 / (y_1-y_0) *2
    center[1] = u[1]+xdiff-1
  end
  if size(S2,2) > 1
    u = floor.(Int, centerOfMass(abs.(S2)))
    println(u)
    y_0 = real(S2[u[1],u[2]-1,u[3]])
    y_1 = real(S2[u[1],u[2]+1,u[3]])
    xdiff = -y_0 / (y_1-y_0) *2
    center[2] = u[2]+xdiff-1
  end
  if size(S3,3) > 1
    u = floor.(Int, centerOfMass(abs.(S3)))
    y_0 = real(S3[u[1],u[2],u[3]-1])
    y_1 = real(S3[u[1],u[2],u[3]+1])
    xdiff = -y_0 / (y_1-y_0) *2
    center[3] = u[3]+xdiff-1
  end
  return center
end


function centerOfMass{T}(data::Array{T,3})
  data = abs.(data)
  maxData = maximum(data)
  data[ data .< maxData*0.1] = 0

  lx=0
  ly=0
  lz=0
  for i3=1:size(data,3)
    for i2=1:size(data,2)
      for i1=1:size(data,1)
           lx += data[i1,i2,i3]*i1
           ly += data[i1,i2,i3]*i2
           lz += data[i1,i2,i3]*i3
       end
     end
  end
  cm = [lx/sum(data), ly/sum(data),lz/sum(data)]
  for d=1:3
    if size(data,d) == 1
      cm[d] = 1.0
    end
  end
  return cm
end
