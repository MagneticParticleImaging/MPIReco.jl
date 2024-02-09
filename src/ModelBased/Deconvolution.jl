export DCR2_deconvolve_l2, DCR2_deconvolve_aniso_l2

function DCR2_deconvolve_l2(c1_2, c2_2, Kern_x, Kern_y, lambda, GradStrx, GradStry)
  # DCR2_deconvolve_l2: Deconvolves DCR with given Langevin-kernel. Input:
  # c1_2, c2_2: Output of DCR2. Kern_x: Kernel used for x-receive-path.
  # lambda: regularization parameter. GradStrx, Gradstry: Gradient strength
  # necessary for correct orientation. 
  
  M, N = size(c1_2);
  
  n_kern_x = size(Kern_x,2);
  n_kern_y = size(Kern_x,1);
  
  n_pad_x = n_kern_x;
  n_pad_y = n_kern_y;
  
  Kern_x_padded = zeros(n_kern_y +n_pad_y, n_kern_x +n_pad_x)
  Kern_x_padded[1:n_kern_y, 1:n_kern_x] .= Kern_x

  Kern_y_padded = zeros(n_kern_y +n_pad_y, n_kern_x +n_pad_x)
  Kern_y_padded[1:n_kern_y, 1:n_kern_x] .= Kern_y

  Kern_x_padded = padarray(Kern_x, Fill(0, (0, 0), (n_pad_x,n_pad_y)))
  Kern_y_padded = padarray(Kern_y, Fill(0, (0, 0), (n_pad_x,n_pad_y)))
  
  #FFTTemp1 = fft2(Kern_x,n_kern_y +n_pad_y,n_kern_x +n_pad_x);
  #FFTTemp2 = fft2(Kern_x.',n_kern_y +n_pad_y,n_kern_x +n_pad_x);
  
  FFTTemp1 = fft(Kern_x_padded)
  FFTTemp2 = fft(Kern_y_padded)
    
   
  L1 = vec(FFTTemp1);
  L2 = vec(FFTTemp2);
  
  Block1 = spdiagm(L1) #spdiags(L1,0,length(L1),length(L1));
  Block2 = spdiagm(L2) #spdiags(L2,0,length(L2),length(L2));

  c1_padded = zeros( n_kern_y +n_pad_y, n_kern_x +n_pad_x)
  c2_padded = zeros( n_kern_y +n_pad_y, n_kern_x +n_pad_x)
  c1_padded[1:size(c1_2,1), 1:size(c1_2,2)] .= sign(GradStrx)*real.(c1_2)
  c2_padded[1:size(c2_2,1), 1:size(c2_2,2)] .= sign(GradStry)*real.(c2_2)

  @info size(c1_padded)

  #c1_padded = padarray(padarray(sign(GradStrx)*real.(c1_2), Pad(:symmetric, (0, 0), (21,21))),Fill(0, (0, 0), (42,42)))
  #c2_padded = padarray(padarray(sign(GradStrx)*real.(c2_2), Pad(:symmetric, (0, 0), (21,21))),Fill(0, (0, 0), (42,42)))

  #c1_padded = padarray(sign(GradStrx)*real.(c1_2), Pad(:symmetric, (0, 0), (63,63)))
  #c2_padded = padarray(sign(GradStrx)*real.(c2_2), Pad(:symmetric, (0, 0), (63,63)))

  #c1_padded = padarray(sign(GradStrx)*real.(c1_2), Pad(:reflect, (0, 0), (n_kern_x+n_pad_x-M,n_kern_y+n_pad_y-N)))
  #c2_padded = padarray(sign(GradStrx)*real.(c2_2), Pad(:reflect, (0, 0), (n_kern_x+n_pad_x-M,n_kern_y+n_pad_y-N)))


  c1_padded = padarray(sign(GradStrx)*real.(c1_2), Fill(0, (0, 0), (n_kern_x+n_pad_x-M,n_kern_y+n_pad_y-N)))
  c2_padded = padarray(sign(GradStrx)*real.(c2_2), Fill(0, (0, 0), (n_kern_x+n_pad_x-M,n_kern_y+n_pad_y-N)))


  @info size(c1_2)
  @info size(c1_padded)

  b1 = fft(c1_padded) #fft2(sign(GradStrx)*real(c1_2),n_kern_y +n_pad_y, n_kern_x +n_pad_x);
  b2 = fft(c2_padded) #fft2(sign(GradStry)*real(c2_2),n_kern_y +n_pad_y, n_kern_x +n_pad_x);
  b = cat(vec(b1), vec(b2), dims=1)

  A = cat(Block1, Block2, dims=1)
  
  b3 = A'*b;
  A2 = A'*A;
  
  A3 = A2 + lambda *spdiagm(ones(size(A2,2))) # speye(size(A2,2));
  
  s = A3\b3;
  
  s2 = ifft(reshape(s,n_kern_y+n_pad_y,n_kern_x+n_pad_x))*4*(n_kern_y+n_pad_y)*(n_kern_x+n_pad_x);

  s2[real.(s2).< 0] .= 0;
  s2 = real.(s2);
  
  
  s2 = s2[end-M+1:end,end-N+1:end];
 
  
  if GradStrx<0
     s2=reverse(s2, dims=2); 
  end
  if GradStry<0
     s2=reverse(s2, dims=1);
  end

  return s2
end



function ConvMatRow(Kern, OL, N_)
  # ConvMatRow creates a row of a convolution matrix. Kern: Matrix containing the kernel. 
  # OL: Position of the upper left corner in the image, where the kernel shall be placed
  # in the image. N_: Size of the image that shall be convolved with the kernel. 
  
  
  N, M = size(Kern);
  
  Kern = reverse(Kern, dims=(1,2)); # due to convolution
  
  Indices = reshape(OL[1] : (OL[1]+N-1),:,1) * ones(Int,1,M) + 
          ones(Int,N,1) * reshape(((OL[2]-1):(OL[2]+M-2)).*N_[1],1,:);
  
  row = zeros(1, prod(N_));
  row[Indices] = vec(Kern);
  return row
end
  
  

function DCR2_deconvolve_aniso_l2(c1_2, c2_2, kernx, kerny, lambda, GradStrx, GradStry, DirIt="direct")
  # DCR2_deconvolve_aniso_l2(c1_2, c2_2, Kern_x, lambda, GradStrx, GradStry, DirIt)
  # Performs a deconvolution of the DCR output with spatially varying kernel.
  # Input: c1_2, c2_2: Output of the DCR-Function. kernels: 3D-Matrix
  # containing the varying kernels. lambda: Regularization parameter.
  # GradStrx, GradStry: Gradient strength, necessary for correct orientation.
  # DirIt: String containing either 'direct' or 'iterative' for choice of the
  # method. 
  
  n_kern_x = size(kernx,2);
  n_kern_y = size(kernx,1);
  
  n_cpre, m_cpre = size(c1_2);
  #c1_pad = padarray(real(c1_2),[(n_kern_y-1)/2, (n_kern_x-1)/2],0,'both'); % actually just needed to compute the size
  size_padded_c1 = [n_cpre+((n_kern_y-1)), m_cpre+((n_kern_y-1))];
  
  
  M_x = zeros(length( c1_2 ), prod(size_padded_c1));
  M_y = zeros(length( c1_2 ), prod(size_padded_c1));
  
  for i = 1:n_cpre
      for j = 1:m_cpre
          #pos = sub2ind([size(c1_2,1), size(c2_2,2)]  , i, j);
          pos = i + (j-1)*size(c1_2,1)
          M_x[i + (j-1)*n_cpre, :] = ConvMatRow(kernx[:,:, pos], [i, j], size_padded_c1);
          M_y[i + (j-1)*n_cpre, :] = ConvMatRow(kerny[:,:, pos], [i, j], size_padded_c1);
      end
  end
  
  
  M = [M_x; M_y];
  b = [real.(vec(c1_2)); real.(vec(c2_2))];
  
  if DirIt == "direct"
  
    A2 = M'*M;
    A2 = A2 + lambda * spdiagm(ones(size(M,2)));
    b3 = M'*b;
    
    
    test = A2\b3;
    test[test.< 0] .= 0;
    tmp = reshape(test, (n_kern_y-1)+n_cpre, (n_kern_x-1)+m_cpre);
  else
    error("No method given")   
  end

  if GradStrx<0
    tmp=reverse(tmp, dims=2); 
 end
 if GradStry<0
    tmp=reverse(tmp, dims=1);
 end
  
  reko_aniso = tmp[(n_kern_y-1)÷2+1:end-(n_kern_y-1)÷2, (n_kern_x-1)÷2+1:end-(n_kern_x-1)÷2];
  return reko_aniso
end

  
