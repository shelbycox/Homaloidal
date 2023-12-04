loadPackage("Graphs")
loadPackage("Matroids")
loadPackage("Cremona")

makeMonomial = (x, v) -> product (for i from 0 to length(entries v)-1 list x_i^(v_i))
matroidPolynomial = (M, x) -> {
    B := basisIndicatorMatrix(M);
    T := for i from 0 to numcols(B)-1 list makeMonomial(x, B_i);
    return sum T
}

E = {{0,1}, {1,2}, {2,3}, {0,3}, {0,4}, {1,4}, {2,4}, {3,4}};
G = Graphs$graph E;
edges G
M = matroid G;
f = matroidPolynomial(M, QQ[x_0..x_7]);
gradf = diff(vars ring f, f);
phi = map(ring f, ring f, gradf);
Phi = rationalMap(phi);
degree Phi
isDominant(phi)
isBirational(phi)