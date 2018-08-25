using Pkg
println("Installing LinearSolver")
Pkg.add(PackageSpec(url="https://github.com/tknopp/LinearSolver.jl.git", rev="julia-0.7"))
println("Installing MPIFiles")
Pkg.add(PackageSpec(url="https://github.com/MagneticParticleImaging/MPIFiles.jl.git", rev="julia-0.7"))
