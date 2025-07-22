# # Imaging Operators
using MPIReco #hide
# The system-matrix based image reconstruction algorithms provided by MPIReco mainly focus on inverse problems of the form:
# ```math
# \begin{equation}
#   \underset{\mathbf{c}}{argmin} \frac{1}{2}\vert\vert \mathbf{S}\mathbf{c}-\mathbf{u} \vert\vert_2^2 + + \mathbf{R(x)}
# \end{equation}
# ```
# where $\mathbf{S}$ is a system matrix, $\mathbf{u}$ is the measurement vector, and $\mathbf{R(x)}$ is an (optional) regularization term.
# In this explanation we will take a closer look at the requirements on $\mathbf{S}$, which is especially relevant if one wants to implement a hybrid- or model based operator for the inverse problem.
S = randn(256, 256)
c = randn(256)
u = S*c;

# ## Linear Algebra
# Different solvers require different interaction with the system matrix $\mathbf{S}$.
# Kaczmarz itself uses a dot product with the rows of $\mathbf{S}$ and the current approximate solution $\mathbf{c}$. Since there is no predefined function for such an operation, RegularizedLeastSquares.jl implemented its own
using MPIReco.RegularizedLeastSquares
row = 1
isapprox(RegularizedLeastSquares.dot_with_matrix_row(S, c, row), sum(S[row, :] .* c))
# Custom operators should implement efficient variants of this function to ensure good performance with the Kaczmarz solver.
# Since Julia is a column-major order language, this row-based access pattern is quite inefficient for dense arrays. A workaround is to transpose the matrix then pass it to a Kaczmarz solver.
S_efficient = transpose(collect(transpose(S)))
typeof(S_efficient)
# RegularizedLeastSquares provides an efficient implmentation for this operator.

# Other solvers such as FISTA or CGNR require the normal operator of $\mathbf{S}$, either because they require the gradient of the least squares norm or because the work on an adapted optimization problem.
# The normal operator $\mathbf{S^*S}$ is composed of a matrix-vector product of $\mathbf{S}$ followed by matrix-vector product of the adjoint of $\mathbf{S}$.
# An efficient matrix-vector product is provided in Julia by the LinearAlgebra standard library:
using LinearAlgebra
mul!(u, S, c)
isapprox(u, S * c)
# This inplace variant of a matrix-vector product is used within in RegularizedLeastSquares and also needs to be provided for custom operators.

# ## Matrix-Free
# LinearAlgebra also provides a function to compute the adjoint of a matrix-vector product:
mul!(c, adjoint(S), u)
isapprox(c, adjoint(S) * u)
# This function creates a lazy adjoint, i.e. it does not create a new array. Instead it only changes the way the size and indexing into the underyling array works.
S_adj = adjoint(S)
# A related concepts are matrix-free operators. These are operators which behave like a matrix in a(n adjoint) matrix-vector product, but don't have an underlying dense matrix.

# There are several different Julia packages with which one can create such matrix-free operators such as LinearOperators.jl or LinearMaps.jl.
# MPIReco uses the [LinearOperatorCollection.jl](https://github.com/JuliaImageRecon/LinearOperatorCollection.jl) package to create matrix-free operators for and with $\mathbf{S}$.
using MPIReco.LinearOperatorCollection
weights = rand(256)
wop = WeightingOp(weights)
size(wop)
# For example, LinearOperatorCollection provides a weighting operator. This operator acts as if it is a diagonal matrix, but only requires the elements of the diagonal.
# We can also do matrix-free matrix product to create $\mathbf{WS}$ without calculating the matrix-matrix product:
WS = ProdOp(wop, S)
isapprox(WS * c, weights .* S * c)
# Another relevant matrix-free operator is the normal operator:
SHS = normalOperator(WS)
# This operator not only computes its matrix-vector product in a lazy fashion, it also optimizes the resulting operator:
# Without an optimization a matrix-free product would apply the following operator each iteration:
# ```math
# \begin{equation}
#   (\mathbf{WS})^*\mathbf{WS} = \mathbf{S}^*\mathbf{W}^*\mathbf{W}\mathbf{S}
# \end{equation}
# ```
# This is not efficient and instead the normal operator can be optimized by initially computing the weights:
# ```math
# \begin{equation}
#   \tilde{\mathbf{W}} = \mathbf{W}^*\mathbf{W}
# \end{equation}
# ```
# and then applying the following each iteration:
# ```math
# \begin{equation}
#   \mathbf{S}^*\tilde{\mathbf{W}}\mathbf{S}
# \end{equation}
# ```

# And lastly, the efficient multi-patch image reconstruction is based on a matrix-free multi-patch operator.
# This operator contains system matrices per patch as well as indexing metadata and is able to act like a large
# dense matrix in a matrix-vector product and can be combined with weighting and the normal operator.

# ## GPU Acceleration
# GPU acceleration is achieved by adapting $\mathbf{S}$ and $\mathbf{u}$ to be GPU-compatible data types. In the case of 
# dense arrays, this means the arrays are GPU arrays.

# Matrix-free operators on the other hand don't necessarily need to be on the GPU. In the case that a matrix-free operator is only
# a function call, this function only needs to be GPU compatible. In the case of the weighting operator or the matrix-free operator,
# we need to move the internally used dense arrays to the GPU and can then use them on the GPU.

# To make more complex operators work on the GPU one can either try to formulate a method on GPUArrays which only uses broadcasts.
# Or alternatively implement a custom GPU kernel with either a specific GPU package such as CUDA.jl or AMDGPU.jl or a vendor-agnostic
# kernel using KernelAbstractions.jl.

# A related Julia package is also the [Adapt.jl](https://github.com/JuliaGPU/Adapt.jl), which offers the `adapt` method. This essentially acts like `convert(T, a)` without the restriction
# that the result must be of type T. The use case for GPU computing here is to provide the GPU array type T and then return a GPU compatibe version of `a`, however that might look like.