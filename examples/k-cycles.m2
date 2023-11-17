loadPackage("Graphs")
loadPackage("Matroids")
loadPackage("Cremona")
loadPackage("NumericalAlgebraicGeometry")

makeMonomial = (x, v) -> product (for i from 0 to length(entries v)-1 list x_i^(v_i))
matroidPolynomial = (M, x) -> {
    B := basisIndicatorMatrix(M);
    T := for i from 0 to numcols(B)-1 list makeMonomial(x, B_i);
    return sum T
}

k = 4 -- very slow for k >= 6
G = Graphs$cycleGraph k
edges G
M = matroid G
f = matroidPolynomial(M, CC[x_0..x_(k-1)])
gradf = diff(vars ring f, f)
phi = map(ring f, ring f, gradf)
Phi = rationalMap(phi)
degree Phi

R = CC[x_0..x_(k-1), y_0..y_(k-1)]
gradf = sub(gradf, ring f)
mat = gradf || (vars R)_{k..(2*k-1)}
mat2 = gradf || matrix {{1,0,1,1}}
I = minors(2, mat2)
use ring f
J = I + ideal (x_0 - 1)
mingens I
dim J
sols = solveSystem J_*

eval = map(ring f, ring f, coordinates sols_3)
eval(gradf)
