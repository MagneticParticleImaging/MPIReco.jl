function download_(filenameServer, filenameLocal)
  if !isfile(filenameLocal)
    @info "download $(filenameLocal)..."
    HTTP.open("GET", "http://media.tuhh.de/ibi/openMPIData/data/"*filenameServer) do http
      open(filenameLocal, "w") do file
        write(file, http)
      end
    end
  end
end

mkpath("./data/")

# Download data
download_("calibrations/17.mdf", "./data/SF.mdf")
download_("measurements/backgroundDrift/1.mdf", "./data/meas.mdf")
