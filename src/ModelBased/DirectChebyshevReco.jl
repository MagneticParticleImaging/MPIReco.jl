export DCR2

function chebyshevU(n,x)
  return sin((n+1)acos(x))/(sqrt(1-x^2))
end

function chebyshevT(n,x)
  return cos(n*acos(x))
end

function Vn(x,n)
  if n!=0
    if abs(x) < 1
      chebyshevU(abs(n)-1,x)*(sqrt(1-x^2))/abs(n)
    else
      0.0
    end
  else
    if x >= 1
     π/2
    elseif x<= 1
     -π/2
    else
      π/2 - acos(x)
    end
  end
end

function freqIndexToChebyshevOrder(K, Nb)
  D = length(K)
  maxK = maximum.(K)

  λ = [ [((2*Nb-1)*k/(2*Nb^2-2*Nb+1)) for k=1:maxK[d]] for d=1:D]
  M = [ cat(-(1:maxK[d])+round.(Int, λ[d])*Nb, (1:maxK[d])-round.(Int, λ[d])*(Nb-1),dims=2) for d=1:D]

  return M, λ
end

function DCR2(f, freqs, N, Nb, fD, excitation="cos", basistype="UU")
  # Two-dimensional direct Chebyshev reconstruction. Input: Vectors f1, f2 with frequency components of
  # voltage signal corresponding to x- and y-receive coil. K12: Vector of
  # frequency components. n_x,n_y: Number of pixels at which the DF-FOV shall
  # be evaluated. Nb: Frequency divider. Excitation: Used excitation, 'sin'
  # or 'cos'. Type: 'UU' or 'UT'. Uses either Chebyshev polynomials (CPs) of second
  # kind only (UU) or a mixed tensor product of CPs of first and second kind
  # (UT). This has influence on the kernels that have to be used for
  # deconvolution. 
  
  D = size(f,2)
  K = [getindex.(filter( d_->d_[2]==d,  freqs), 1).-1 for d=1:D]
  
  M, λ = freqIndexToChebyshevOrder(K, Nb)

  ###
  M_ = [ (abs(M[1][l,1]), abs(M[1][l,2])) for l=1:size(M[1],1)]
  Mu = unique(M_)
  kl = [ findlast(d->d == Mu[k],M_) for k=1:length(Mu)]

  Q = zeros(Int, size(M[1]))
  Q[kl,:] .= M[1][kl,:]
  M[1] = M[2] = Q



  c = zeros(ComplexF64, N..., D);
  
  # different basis for different excitation
  if excitation == "sin"
      excterm = (-1im);
      factor = 1;
  elseif excitation == "cos"
       excterm = (-1);
       factor = 1/1im;
  else
      error("No valid excitation")
  end
              
  eps_ = 0.01
  x_ = [range(-1+eps_,1-eps_,length=N[d]) for d=1:D]
  
  # Use CPs of 2nd kind only
  if basistype == "UU"
    for d=1:D
        for k in K[d]
          if all(M[d][k,:].!=0) 
            for x in CartesianIndices(ntuple(d->(1:N[d]), D))
              temp = 4/2/pi/fD/(-pi^2) *factor/k/excterm.^(round(Int, λ[d][k])+1) * abs(prod(M[d][k,:])) 
              for t=1:D
                temp *= chebyshevU(abs(M[d][k,t])-1, x_[t][x[D-t+1]]) 
              end
              c[x, d] += f[k+1,d]*temp
            end
          end
        end
      end


    # Use tensorproduct of CPs of 1st and 2nd kind
    elseif basistype == "UT"
      for d=1:D
        for k in K[d]
            if all(M[d][k,:].!=0) 
              for x in CartesianIndices(ntuple(d->(1:N[d]), D))
                temp = 4/2/pi/fD/(-pi^2) *factor/k/excterm.^(round(Int, λ[d][k])+1) * abs(prod(M[d][k,:])) 
                for t=1:D
                  if t==d
                    temp *= chebyshevU(abs(M[d][k,t])-1, x_[t][x[D-t+1]])
                  else
                    temp *= chebyshevT(abs(M[d][k,t]), x_[t][x[D-t+1]])
                  end
                end
                c[x, d] += f[k+1,d]*temp
            end
          end
        end
      end
    end

  return c 
end



function chebyshevForwardModel(K_, N, Nb, excitation="cos")

  D = length(N)
  K = 1:K_
  
  M, λ = freqIndexToChebyshevOrder(repeat([K],D), Nb)

  # different basis for different excitation
  if excitation == "sin"
      excterm = (-1im);
      factor = 1;
  elseif excitation == "cos"
        excterm = (-1);
        factor = 1/1im;
  else
      error("No valid excitation")
  end
              
  eps_ = 0.01
  x_ = [range(-1+eps_,1-eps_,length=N[d]) for d=1:D]
    
  S = zeros(ComplexF64, N..., K, D);
  
  for d=1:D
    for k in 2:K_ # K
      for x in CartesianIndices(ntuple(d->(1:N[d]), D))
        temp = 2/pi * factor / (k-1) / excterm.^(round(Int, λ[d][k-1])+1) 
        for t=1:D
          temp *= Vn(x_[t][x[D-t+1]], M[d][k-1,t]) 
        end
        S[x, k, d] = temp
      end
    end
  end
  
  return S
end

function convCheb(S, Kernx, Kerny)
  Sout = similar(S)
  for k = 1:size(S,3)
    Sout[:,:,k,1] =    imfilter(real(S[:,:,k,1]), Kernx) + 
                    im*imfilter(imag(S[:,:,k,1]), Kernx)

  end
  for k = 1:size(S,3)
    Sout[:,:,k,2] =    imfilter(real(S[:,:,k,2]), Kerny) + 
                    im*imfilter(imag(S[:,:,k,2]), Kerny)

  end
  return Sout
end