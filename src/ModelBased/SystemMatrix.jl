export calcSM

function direction(x)
  if norm(x) == 0
    return (1.0,0.0,0.0)
  else
    return x ./ norm(x)
  end
end

function calcSM(params; chebyshev=false)

  amplitude = params[:amplitude]
  gradient = params[:gradient] 
  grid = params[:grid]

  dividers = params[:dividers]
  tLengthSM = lcm(dividers...);             # length of time vector
  tMaxSM = tLengthSM / params[:samplingRate]; # maximum evaluation time in seconds
  
  tSM = range(0, step=1/params[:samplingRate], length=tLengthSM);
  
  freq = ntuple( d -> dividers[d] > 1 ? params[:samplingRate] ./ dividers[d] : 0.0, 3)
  ampl = ntuple( d -> dividers[d] > 1 ? amplitude[d] : 0.0, 3)

  BSM = (t, offset) -> MNPDynamics.SVector{3,Float32}(ampl[1]*cos(2*pi*freq[1]*t)+offset[1], 
                                          ampl[2]*cos(2*pi*freq[2]*t)+offset[2], 
                                          ampl[3]*cos(2*pi*freq[3]*t)+offset[3] )
  nOffsets = shape(grid)
  factor = params[:kAnis]
  kAnisγ = params[:kAnisγ]
  anisotropyAxis = get(params, :anisotropyAxis, nothing)

  offsets = vec([ gradient.*Tuple(x)  for x in grid ])

  dirs = vec([ (direction(Tuple(x)).*(norm.(x).^kAnisγ))  for x in grid ])
  maxDir = maximum(norm.(dirs))
  
  if anisotropyAxis == nothing
    anisotropyAxis = vec([ (factor./maxDir).*(direction(Tuple(x)).*(norm.(x).^kAnisγ))  for x in grid ])
  else
    anisotropyAxis = vec([ Tuple(factor.*anisotropyAxis)  for x in grid ])
  end
  params_ = copy(params)
  params_[:kAnis] = anisotropyAxis

  if chebyshev
    sm = calcSMReducedEq(params_, BSM, tSM, offsets)
    sm .*= reshape(2*pi*im.*(0:(size(sm,3)-1)),1,1,:)
  else
    sm = simulationMNPMultiParams(BSM, tSM, offsets; params_...)
    sm = rfft(reshape(permutedims(sm,(3,1,2)), nOffsets[1:2]..., :, 3), 3);
    sm .*= reshape(2*pi*im.*(0:(size(sm,3)-1)),1,1,:) / ((size(sm,3)-1)*2)
  end
  return sm
end



function calcSMReducedEq(params, BSM, tSM, offsets)
  
  amplitude = params[:amplitude]
  gradient = params[:gradient] 
  grid = params[:grid]

  #drive-field has cosine-exciation, otherwise sine-exciation
  isCosine = true;

  Xcord, Ycord, Zcord = ([ x[d] for x in grid ] for d=1:3)
  
  kB = 1.38064852e-23
  gamGyro = 1.75*10.0^11
  μ₀ = 4*π*1e-7 #vacuum permeability
  MS = 474000.0;
  temp = 293.0;
  
  fB = params[:samplingRate]
  dividers = params[:dividers]
  tLengthSM = lcm(dividers...);             # length of time vector
  TD = tLengthSM / params[:samplingRate]; # maximum evaluation time in seconds

  Nb = round(Int,lcm(params[:dividers]...) / params[:dividers][2]) #17; # frequency divider
  # frequency index to calculated
  maxFreqIdx = round(Int, TD * params[:samplingRate])÷2+1
  K = 0:(maxFreqIdx-1);
  # find optimal lambda_k values see "10.1088/1361-6560/ac4c2e"
  lambdak = round.(Int, (2*Nb-1)*K/(2*Nb^2-2*Nb+1));
  # relate the lambda_k to the mixing order
  nk = -K+lambdak*Nb;
  mk = K-lambdak*(Nb-1);
  
  # please do not change the sign of the selection field, otherwise the phase
  # must be manipulated

  GradStrx = gradient[1];
  GradStry = gradient[2];
  GradStrz = gradient[3];
  
  # amplitudes can be manipulated.
  Ax = amplitude[1]; # amplitude | x-direction drive field
  Ay = amplitude[2]; # amplitude | y-direction drive field
  Az = amplitude[3]; # amplitude | z-direction
  
  # To be able to generate the time derivative of the mean magentic moment
  omegak = 2*pi*K'/TD*(1im);
  
  
  # phase term to obtain a sinusoidal drive field excitation
  phase_ = (1im).^(lambdak)/pi^2;
  if isCosine
      # phase manipluation to obatin a cosine excitation
      phase_ = phase_.*exp.(1im.*pi./2.0*lambdak);
  end
  
  
  order = 120; # order factor for truncation of series expansion -> large values could lead to NaNs or Inf.
  # small order faster calculation, but lower accuracy
  epsilon = 1e-10; #to avoid divisions by values close to zero in the series expansion
  
  
  nEasy_H_S = collect([ (norm(x) > 0 ? x[d]./norm(x) : (d==1)) for x in params[:kAnis], d=1:3]')
  g_Aniso = norm.(params[:kAnis])
  
  
  #Density the for numerical integration with Chebyshev grid 
  nChebGrid = 51; #201;
  
  #Create a Chebyshev grid around each sampling grid point.
  xk = Ax./GradStrx.*cos.(pi.*((1:nChebGrid).-1/2)./nChebGrid);
  yk = Ay./GradStry.*cos.(pi.*((1:nChebGrid).-1/2)./nChebGrid);
  
  # Coordinates that must be evaluated and padding of K_Aniso and the easy axis
  deltaX, deltaY = meshgrid(xk,yk);
  X_cords_gauss_quad = zeros(length(xk),length(yk),length(Xcord));
  Y_cords_gauss_quad = zeros(length(xk),length(yk),length(Xcord));
  Z_cords_gauss_quad = zeros(length(xk),length(yk),length(Xcord));
  g_aniso_quad = zeros(length(xk),length(yk),length(Xcord));
  nEasy_quad = zeros(3,length(xk),length(yk),length(Xcord));
  
  
  for i=1:length(Xcord)
      Xq = Xcord[i].-deltaX;
      Yq = Ycord[i].-deltaY;
      Zq = Zcord[i].*ones(size(deltaX));
      X_cords_gauss_quad[:,:,i] = Xq;
      Y_cords_gauss_quad[:,:,i] = Yq;
      Z_cords_gauss_quad[:,:,i] = Zq;
      g_aniso_quad[:,:,i] .= g_Aniso[i];
      nEasy_quad[:,:,:,i] = repeat(nEasy_H_S[:,i],1,size(deltaX,1),size(deltaX,2));
  end

  
  # The magentic field to be evaluated
  H = hcat(GradStrx*X_cords_gauss_quad[:],GradStry*Y_cords_gauss_quad[:],GradStrz*Z_cords_gauss_quad[:]);

  #calculated magnetic moment of the given field H
  Mag = MNPDynamics.eqAnisoMeanMagFullVec(H, params[:DCore], MS, temp, vec(g_aniso_quad), 
                                          reshape(nEasy_quad,3,:), order, epsilon);
    
  Mag = reshape(Mag,size(X_cords_gauss_quad)..., 3);

  
  # Perform the convolution using Chebyshev nodes as DCT-II
  Mag_ = dct(Mag, 1:2);
  
  # adjust the weighting of Matlab DCT-II
  weightsY = ones(size(Mag_,1),1)*(pi/sqrt(nChebGrid));
  weightsY[2:end] = weightsY[2:end]./sqrt(2);
  weightsX = weightsY';
  
  
  Mag_ = (Mag_.*weightsY).*weightsX;
  Mag_ = reshape(Mag_,size(Mag_,1),size(Mag_,2),length(grid),3);
  
  S = zeros(ComplexF64, size(Mag_,3),length(K), 3);
  
  
  # Using the DCT coefficients directly as mixing order. The result
  # corresponds to the convolution between weighted Chebyshev polynomials and \varepsilon
  for ij=1:length(K)
      S[:,ij,1] = phase_[ij]*vec(Mag_[abs(mk[ij])+1,abs(nk[ij])+1,:,1]);
      S[:,ij,2] = phase_[ij]*vec(Mag_[abs(mk[ij])+1,abs(nk[ij])+1,:,2]);
      S[:,ij,3] = phase_[ij]*vec(Mag_[abs(mk[ij])+1,abs(nk[ij])+1,:,3]);
  end
  
  
  return reshape(S, shape(grid)[1], shape(grid)[2], length(K), 3)
end