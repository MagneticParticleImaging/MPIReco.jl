using HTTP

# Doanload data
if !isdir("data")
  @info "download data.zip"
  HTTP.open("GET", "http://media.tuhh.de/ibi/MotionCompensation/data.zip") do http
    open("data.zip", "w") do file
        write(file, http)
    end
  end
  @info "extracting data.zip"
  run(`unzip -oq data.zip`)
  rm("data.zip")
end

if !useCompressedMatrices && !isfile("data/SF1Large.mdf")
    @info "download data.zip"
    HTTP.open("GET", "http://media.tuhh.de/ibi/MotionCompensation/dataLarge.zip") do http
      open("dataLarge.zip", "w") do file
          write(file, http)
      end
    end
    @info "extracting dataLarge.zip"
    run(`unzip -oq dataLarge.zip`)
    rm("dataLarge.zip")
end
