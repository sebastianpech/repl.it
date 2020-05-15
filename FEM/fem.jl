using Pkg
pkg"activate ."
using ForwardDiff
using Plots
using LinearAlgebra

## # Element geometry
# 2D linear element with the following coodinates

X = [[ -1,  1 ], [  0.5,  1 ], [  1,  1 ], 
     [ -1,  0 ],             [  1,  0.5 ], 
     [ -1, -1 ], [  0, -1 ], [  1, -1 ]]

## # Shape function definition
# First one needs the geometry in natural coordinates:

Xn = [[ -1,  1 ], [  0,  1 ], [  1,  1 ], 
      [ -1,  0 ],             [  1,  0 ], 
      [ -1, -1 ], [  0, -1 ], [  1, -1 ]]

# The shape functions can be computed by fitting a polynomal function,
# such that it is only 1 at the node the shape function is generated on
# and 0 on all other nodes

npoly(ncoeff, x...) = npoly(Val(ncoeff), x...)
npoly(ncoeff, x::Vector) = npoly(Val(ncoeff), x...)
npoly(::Val{2}, x) = [x^0 x]
npoly(::Val{3}, x) = [x^0 x x^2]
npoly(::Val{4}, x, y) = [x^0, x, y, x*y]
npoly(::Val{8}, x, y) = [x^0, x, y, x*y, x^2, y^2, x*y^2, x^2*y]

function generate_shape_functions(Xn,j)
    ncoeff = length(Xn)
    # Compute coefficient matrix for all coordinates
    A = vcat(transpose.(npoly.(ncoeff,Xn))...)
    # The function is only allowed to be unequal to zero at node j
    F = zeros(ncoeff)
    F[j] = 1.0
    # Solve for coefficients
    a = A\F
    (ξ, η)->( npoly(ncoeff, ξ, η)⋅a )
end

N = [ generate_shape_functions(Xn, i) for i in 1:length(Xn) ]

## # Jacobi matrix
# Compute the jacobi matrix by first defining the real coordinates in term
# of the natural coodinates

# Definition of isoparametric elements
x(xn) = sum( N[i](xn...)*X[i] for i in 1:length(X) )

# Definition of jacobi matrix
J = xn -> transpose(ForwardDiff.jacobian(x, xn,))

surface(-1:0.1:1, -1:0.1:1, (ξ, η)->(det(J([ ξ,η ]))))
