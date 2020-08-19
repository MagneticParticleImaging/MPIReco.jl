using HTTP

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
download_("calibrations/12.mdf", "./data/SF1Small.mdf")
download_("calibrations/13.mdf", "./data/SF2Small.mdf")
download_("calibrations/14.mdf", "./data/SF3Small.mdf")
download_("calibrations/15.mdf", "./data/SF4Small.mdf")

if useFastData
  download_("measurements/rotationPhantom/2.mdf", "./data/measFast.mdf")
else
  download_("measurements/rotationPhantom/1.mdf", "./data/measSlow.mdf")
end
download_("measurements/rotationPhantom/3.mdf", "./data/measBG.mdf")
download_("measurements/rotationPhantom/4.mdf", "./data/measStatic.mdf")

if !useCompressedMatrices
  download_("calibrations/8.mdf", "./data/SF1Large.mdf")
  download_("calibrations/9.mdf", "./data/SF2Large.mdf")
  download_("calibrations/10.mdf", "./data/SF3Large.mdf")
  download_("calibrations/11.mdf", "./data/SF4Large.mdf")
end
