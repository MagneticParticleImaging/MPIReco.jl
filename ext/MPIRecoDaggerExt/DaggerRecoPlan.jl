function Base.setproperty!(dplan::DaggerImageReconstruction.DaggerRecoPlan{<:AbstractMPIRecoParameters}, name::Symbol, data::MPIFiles.DMPIFile)
  fetch(Dagger.spawn(dplan._chunk) do chunk
    # If we are on the same worker we can just get the file
    if myid() == MPIFiles.worker(data)
      data = collect(data.file)
    end
    Base.setproperty!(chunk, name, data)
    return true
  end)
end


function Base.setproperty!(dplan::RecoPlan{<:AbstractMPIRecoParameters}, name::Symbol, data::MPIFiles.DMPIFile)
  # If we are on the same worker we can just get the file
  if myid() == MPIFiles.worker(data)
    data = collect(data.file)
  end
  Base.setproperty!(dplan, name, data)
  return data
end
