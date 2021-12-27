

###############

function loss(model, S, c, u, λ, cTrue)
    return Flux.Losses.mse(S*vec(model(c)), u) +  λ*Flux.Losses.mse(vec(model(c)), vec(cTrue)) 
  end
  
  
  ###############
  
  function reconstructionDIP(S, N, u; model = make_model_unet(n), 
                                      opt = ADAM(0.001, (0.7, 0.999)),
                                      cTrue = ones(Float32,N...),
                                      noiseData = vec(rand(Float32,N...).*Float32(0.7)),
                                      iterations = 100,
                                      λ = Float32(0),
                                      iters = nothing)
  
    parameters = params(model)
    data = [(model, S, noiseData, u, λ, device(cTrue))]; # λ, device(cTrue)
    
    @info "training"
    iter = 1:iterations
    #@showprogress 1.0 "Computing..." 
    residual = Float64[]
    error = Float64[]
    ssims = Float64[]
    for i=iter
      Flux.train!(loss, parameters, data, opt)
      currentC = model(noiseData) |> cpu
      push!(residual, loss(data[1]...))
      push!(error, norm(vec(currentC)-vec(cTrue))/norm(vec(cTrue)) )
      push!(ssims, ssim(squeeze(currentC),squeeze(cTrue)) )
      if iters != nothing
        push!(iters, currentC)
      end
      if mod(i,20) == 1
        @info "res=$(residual[end]) err=$(error[end]) ssim=$(ssims[end]) time=$(i/iter[end])"
      end
    end
    cTrained = model(noiseData) |> cpu
  
    return cTrained, residual, error, ssims
  end
  
  
  
  
  function lossPre(model, c, u)
    return Flux.Losses.mae(vec(model(c)), vec(u))
  end
  
  
  ###############
  
  function reconstructionDIPPre(N; model = make_model_unet(n), 
                                      opt = ADAM(0.001, (0.7, 0.999)),
                                      cTrue = ones(Float32,N...),
                                      noiseData = vec(rand(Float32,N...).*Float32(0.7)),
                                      iterations = 100)
  
    parameters = params(model)
    data = [(model, noiseData, cTrue)];
    
    @info "training"
    iter = 1:iterations
    #@showprogress 1.0 "Computing..." 
    residual = Float64[]
  
    for i=iter
      Flux.train!(lossPre, parameters, data, opt)
  
      push!(residual, lossPre(data[1]...))
  
      if mod(i,20) == 1
        @info "res=$(residual[end])  time=$(i/iter[end])"
      end
    end
    cTrained = model(noiseData) |> cpu
  
    return cTrained, residual
  end
  
  
  
  
  