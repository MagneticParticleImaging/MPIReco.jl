export findCenterOfDfFov

type SFTuple
	SF1::Array{Complex{Float64},3}
	SF2::Array{Complex{Float64},3}
    SF3::Array{Complex{Float64},3}
end

@doc "This function searches the middle zero of the vector vec.\n
	For this at first it searches the global extrema of the vector and afterwards the zero between these extrema.\n
	The zero is calculated on subpixels via linear interpolation."->

function searchZero(vec)
    if size(vec)[1]==1
        return 1
    end
    max,maxI=findmax(vec)
    min,minI=findmin(vec)
    if minI<maxI
        x=Int(round(minI+(maxI-minI)/2))
        d=1
    else
        x=Int(round(maxI+(minI-maxI)/2))
        d=-1
    end
	while vec[x]*vec[x+1]>0
		vec[x]>0 ? x=x-d : x=x+d
	end
    # linear interpolation between x and x+1
    a=vec[x+1]-vec[x]
    b=vec[x]-a*x
    return -Float64(b)/a
end

@doc "This function reads the component belonging to the given mixing-factors und chanal from the given MPIfile.\n
In the case of a 2D systemfunction, the third dimension is added." ->

function readSF(bSF,mx::Int,my::Int,mz::Int,recChan::Int,maxFreq::Int)
	freq = mixFactorToFreq(bSF, mx, my, mz)
	freq = clamp(freq,0,maxFreq-1)
	k = freq + (recChan-1)*maxFreq
	SF = getSF(bSF,Int64[k+1],returnasmatrix = false, bgcorrection= true)
    # add third dimension if SF is 2D
    if ndims(SF)==2
        SF=reshape(SF,size(SF)[1],size(SF)[2],1)
    end
    return SF
end

@doc "This function displays the result in th 2D-case"->

function dispResult2D(SF,c::Array,wDFFov,sp)
    if sp==1 figure() end
    subplot(1,2,sp)

    PyPlot.gray()
    imshow(abs(SF[:,:,Int(round(c[3]))])',interpolation="nearest")
    xticks([])
    yticks([])
	x=[c[1]-1,c[1]-1]
    y=[c[2]-4,c[2]+2]
    plot(x,y,color="r")
	x=[c[1]-4,c[1]+2]
	y=[c[2]-1,c[2]-1]
	plot(x,y,color="r",)

	x=[c[1]-1-wDFFov[1],c[1]-1+wDFFov[1],c[1]-1+wDFFov[1],c[1]-1-wDFFov[1],c[1]-1-wDFFov[1]]
	y=[c[2]-1-wDFFov[2],c[2]-1-wDFFov[2],c[2]-1+wDFFov[2],c[2]-1+wDFFov[2],c[2]-1-wDFFov[2]]
	plot(x,y,color="g")
end

@doc "This function displays the result in th 3D-case"->

function dispResult3D(SF,c::Array,wDFFov)
    figure()
    grid=size(SF)[1]

    z=Array(Int,grid^3)
    y1=Array(Int,grid^2)
    for i=1:grid
        z[(i-1)*grid^2+1:i*grid^2]=i
        y1[(i-1)*grid+1:i*grid]=i
    end
    y=repmat(y1,grid)
    x=repmat(1:grid,grid^2)

    d=reshape(abs(SF),grid^3)
	max,=findmax(d)
	p=[x y z].*(d.>0.4*max)
	# delete zeros
	xc=p[:,1][find(p[:,1])]
	yc=p[:,2][find(p[:,2])]
	zc=p[:,3][find(p[:,3])]

    gca(projection="3d")
    scatter(xc,yc,zs=zc,s=1)
    scatter(c[1],c[2],zs=c[3],color="r")

    # green border
    x=[c[1]-1-wDFFov[1],c[1]-1+wDFFov[1],c[1]-1+wDFFov[1],c[1]-1-wDFFov[1],c[1]-1-wDFFov[1],c[1]-1-wDFFov[1]]
	y=[c[2]-1-wDFFov[2],c[2]-1-wDFFov[2],c[2]-1+wDFFov[2],c[2]-1+wDFFov[2],c[2]-1-wDFFov[2],c[2]-1-wDFFov[2]]
    z=[c[3]-1-wDFFov[3],c[3]-1-wDFFov[3],c[3]-1-wDFFov[3],c[3]-1-wDFFov[3],c[3]-1-wDFFov[3],c[3]-1+wDFFov[3]]
	plot3D(x,y,z,color="g")

    x=[c[1]-1+wDFFov[1],c[1]-1+wDFFov[1],c[1]-1-wDFFov[1],c[1]-1-wDFFov[1],c[1]-1+wDFFov[1],c[1]-1+wDFFov[1]]
	y=[c[2]-1+wDFFov[2],c[2]-1-wDFFov[2],c[2]-1-wDFFov[2],c[2]-1+wDFFov[2],c[2]-1+wDFFov[2],c[2]-1+wDFFov[2]]
    z=[c[3]-1+wDFFov[3],c[3]-1+wDFFov[3],c[3]-1+wDFFov[3],c[3]-1+wDFFov[3],c[3]-1+wDFFov[3],c[3]-1-wDFFov[3]]
	plot3D(x,y,z,color="g")

    x=[c[1]-1+wDFFov[1],c[1]-1+wDFFov[1]]
    y=[c[2]-1-wDFFov[2],c[2]-1-wDFFov[2]]
    z=[c[3]-1+wDFFov[3],c[3]-1-wDFFov[3]]
    plot3D(x,y,z,color="g")

    x=[c[1]-1-wDFFov[1],c[1]-1-wDFFov[1]]
    y=[c[2]-1+wDFFov[2],c[2]-1+wDFFov[2]]
    z=[c[3]-1+wDFFov[3],c[3]-1-wDFFov[3]]
    plot3D(x,y,z,color="g")

end

@doc "This function calculates the center of a component-tuple with the initial solution cStart"->

function calcCenter(SF::SFTuple,cStart)
	c=cStart
	# repeat twice to get better result
	# first step selects the slice that is taken in the second step
	for step=1:2
        vec1=real(SF.SF1[:,Int(round(c[2])),Int(round(c[3]))])
        vec2=real(SF.SF2[Int(round(c[1])),:,Int(round(c[3]))])
        vec3=real(SF.SF3[Int(round(c[1])),Int(round(c[2])),:])
        c[1]=searchZero(vec1)
        c[2]=searchZero(vec2)
        c[3]=searchZero(vec3)
	end
	return c
end

@doc "helper fot calculation via zeros"->

function calcWithZero!(bSF,nrOfTuples,c0,maxFreq,combine)
    combine ? c=repmat(c0,1,2*nrOfTuples^3) : c=repmat(c0,1,2*nrOfTuples)
    for i=1:2*nrOfTuples
        chan=Int(ceil(i/nrOfTuples))
        combine ? max=nrOfTuples : max=1
        for j=1:max
            for k=1:max
                a2=i-(chan-1)*nrOfTuples
                if combine
                    a1=j
                    a3=k
                else
                    a1=a2
                    a3=a2
                end
                if chan==1
                    SF=SFTuple(readSF(bSF,2*a1,2,0,chan,maxFreq),readSF(bSF,1,2*a2+1,0,chan,maxFreq),readSF(bSF,1,2,2*a3-1,chan,maxFreq))
                else
                    SF=SFTuple(readSF(bSF,2*a1+1,1,0,chan,maxFreq),readSF(bSF,2,2*a2,0,chan,maxFreq),readSF(bSF,2,1,2*a3-1,chan,maxFreq))
                end
                c[:,(i-1)*max^2+(j-1)*max+j]=calcCenter(SF,c0)
            end
            #println("Chan:",chan," Pair:",k1," ",k2," c:",c[:,(i-1)*maxJ+j])
        end
    end
    return c
end

@doc "helper fot calculation via center of mass"->

function calcWithCOM!(bSF,depth,grid,maxFreq)
    c=zeros(Float64,3,2*depth^3)
    for chan=1:2
        chan==1 ? offset=[0,1,2] : offset=[1,0,2]
        for i=1:depth
            for j=1:depth
                for k=1:depth
                    if i%2==1 && j%2==1
                        nrOfMax=1
                    elseif i%2==0 && j%2==0
                         nrOfMax=4
                    else
                        nrOfMax=2
                    end
                    if k%2!=0 && grid[3]>1
                        nrOfMax=nrOfMax*2
                    end
                    SF=readSF(bSF,i+offset[1],j+offset[2],k+offset[3],chan,maxFreq)
                    center,=calcCenterOfMass_AE(abs(SF),nrOfMax,0.3)
                    pos=(chan-1)*depth^3+(i-1)*depth^2+(j-1)*depth+k
                    c[:,pos]=mean(center,2)
                    #println("Chan:",chan," Komp:",i," ",j," ",k," c:",c[:,pos])
                end
            end
        end
    end
    return c
end

@doc "This function calculates the center of the DfFov of the given MPIFile bSF\n
	dispRes: graphic display of result\n
	combine: combine inside the tuples to get more sets for calculation with same depth\n
	useCOM: true: calculate via center of mass, false: calculate via zeros\n
	depth: number of tuples for searching the zeros or depth of the grid for searching the center of mass\n
	sp: choose subplot in case of 2D-plotting of result"->

function findCenterOfDfFov(bSF::MPIFile;dispRes::Bool=false,combine::Bool=false,useCOM=false,depth::Int=3,sp::Int=1)
	maxFreq=length(frequencies(bSF))
    grid=gridSize(bSF)
	dfFovSize=dfFov(bSF)
	fovSize=fov(bSF)
	pixelsDfFov=ceil.(grid.*dfFovSize./fovSize)
    SFres=readSF(bSF,2,3,0,1,maxFreq)

    if !useCOM
        # calculates initial value by center of mass
        center,=calcCenterOfMass_AE(abs.(SFres),4,0.3)
        c0=mean(center,2)
        c=c0#calcWithZero!(bSF,depth,c0,maxFreq,combine)
    else
        c=calcWithCOM!(bSF,depth,grid,maxFreq)
    end

    cmean=mean(c,2)
    cvar=var(c,2)
    println("Result: ",cmean)
    println("Variance: ",cvar)

	if dispRes
		wDFFov=convert(Array{Int},(floor(pixelsDfFov./2)))
        grid[3]==1 ? dispResult2D(SFres,cmean,wDFFov,sp) : dispResult3D(SFres,cmean,wDFFov)
	end
    return cmean
end

function filterCutOff_{T}(imagestack::AbstractArray{T,3}; relThresh=0.3)
  m = repeat(maximum(imagestack,(1,2,3))[:], outer = [prod(size(imagestack)[1:3])])
  m = reshape(m, size(imagestack))
  return (imagestack .> relThresh*m) .* imagestack
end

@doc "Performs a cutoff filter operation with `relThresh` on `data3D`.
Algorithm tries to find a number of local maxima `numOfMaxima` in the data3D
and calcutes their center of mass: centerOfMass [x,y,z] columnwise"->
function calcCenterOfMass_AE{T}(data3D::Array{T,3}, numOfMaxima, relThresh)
  cutOff=filterCutOff_(data3D, relThresh=relThresh)
  # Convert to binary image
  indecies=find(cutOff)
  cutOff[indecies]=1
  # initialize
  centerOfMass=zeros(eltype(data3D),(3,numOfMaxima))
  xS,yS,zS=size(data3D)
  centerdata3D=zeros(eltype(data3D), xS, yS, zS, numOfMaxima)
  # loop for search of numOfMaxima local maxima
  for k=1:numOfMaxima
    components3D=zeros(data3D)
    # get maximum voxel value
    maxValues, indicesMax = findmax(data3D,[1,2,3])
    # segment connecting areas
    comp=label_components(cutOff,[1,2,3])
    compindices=Images.component_indices(comp)
    # find connection area where the maximum value is in
    compBrightes=comp[indicesMax[:]]
    compBrightesIndices = compindices[compBrightes+1]
    components3D[compBrightesIndices[1]]=data3D[compBrightesIndices[1]]
    centerOfMass[:,k]=centerOfMass_AE(components3D, compBrightesIndices[1])
    # set all connecting area of the local maximum to zero to find the next local maximum
    # and not the next maximum voxel within the area of the previous
    centerdata3D[:,:,:,k]=components3D
    data3D[compBrightesIndices[1]]=0.0
  end
  return centerOfMass, centerdata3D
end

@doc "Performs a cutoff filter operation with `relThresh` on every frame of `data4D`.
Algorithm tries to find a number of local maxima `numOfMaxima` in the each frame of data4D
and calcutes their center of mass"->
function calcCenterOfMass_AE{T}(data4D::Array{T,4}, numOfMaxima, relThresh)
  numFrames=size(data4D,4)
  centers=zeros(eltype(data4D),3,numOfMaxima, numFrames)
  for k=1:numFrames
    cs,unused=calcCenterOfMass_AE(data4D[:,:,:,k], numOfMaxima, relThresh)
    centers[:,:,k]=cs
  end
 return centers
end

function centerOfMass_AE(data::Matrix)
  lx=0
  ly=0
  for i2=1:size(data,2)
      for i1=1:size(data,1)
           lx += abs(data[i1,i2])*i1
           ly += abs(data[i1,i2])*i2
       end
  end
  return [lx/sum(abs(data)), ly/sum(abs(data))]
end

function centerOfMass_AE{T}(data::Array{T,3})
  data = abs(data)
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
  return [lx/sum(data), ly/sum(data),lz/sum(data)]
end

function centerOfMass_AE{T}(data::Array{T,3}, linIndices)
  cm=zeros(T,3)
  for i in linIndices
      iX,iY,iZ=ind2sub(data,i)
      cm+=data[i]*[iX, iY, iZ]
  end
  cm=cm./sum(data[linIndices])
  return cm
end
