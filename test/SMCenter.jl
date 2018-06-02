using MPIReco
using PyPlot

SFs = ["SF_MP01", "SF_MP02", "SF_MP03", "SF_MP04"]

figure(1)
clf()

for (l,SF) in enumerate(SFs)
  subplot(2,2,l)
  bSF = MPIFile(SF)
  a = findCenterOfDfFov(bSF)
  println(a)
  S = MPIReco.readSF(bSF,2,0,0,2);
  imshow(abs.(S)[:,:,1]);
  plot([a[2]-1],[a[1]-1],"gx",lw=4)
end

figure(2)
clf()

for (l,SF) in enumerate(SFs)
  subplot(2,2,l)
  bSF = MPIFile(SF)
  a = findCenterOfDfFov(bSF)
  println(a)
  S = MPIReco.readSF(bSF,0,2,0,1);
  imshow(abs.(S)[:,:,1]);
  plot([a[2]-1],[a[1]-1],"gx",lw=4)
end
