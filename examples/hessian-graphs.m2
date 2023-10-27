loadPackage "Graphs"
loadPackage "Matroids"

makeMonomial = (x, v) -> product (for i from 0 to length(entries v)-1 list x_i^(v_i))
matroidPolynomial = (M, x) -> {
    B := basisIndicatorMatrix(M);
    T := for i from 0 to numcols(B)-1 list makeMonomial(x, B_i);
    return sum T
}

G = Graphs$graph {{0,1}, {1,2}, {2,3}, {0,3}}
M = matroid G
R = QQ[x_(0,1), x_(0,3), x_(1,2), x_(2,3)]
f = matroidPolynomial(M, R_*)

Jac = diff(vars R, f)
H = diff(transpose(vars R), Jac)
h = det H
factor h -- the hessian of a cycle graph is non-vanishing :(