loadPackage("Graphs")
loadPackage("Matroids")
loadPackage("Cremona")

makeMonomial = (x, v) -> product (for i from 0 to length(entries v)-1 list x_i^(v_i))
matroidPolynomial = (M, x) -> {
    B := basisIndicatorMatrix(M);
    T := for i from 0 to numcols(B)-1 list makeMonomial(x, B_i);
    return sum T
}

-- 4-cycle summed with a traingle along an edge
E = {{0,1}, {1,2}, {2,3}, {0,3}, {0,4}, {1,4}};
G = Graphs$graph E;
edges G
M = matroid G;
f = matroidPolynomial(M, QQ[x_0..x_(#E-1)])
gradf = diff(vars ring f, f)
phi = map(ring f, ring f, gradf);
Phi = rationalMap(phi);
isDominant(phi)
isBirational(phi)
degree Phi

-- 5-cycle summed with an edge along a vertex
E = {{0,1}, {1,2}, {2,3}, {3,4}, {0,4}, {0,5}};
G = Graphs$graph E;
edges G
M = matroid G;
f = matroidPolynomial(M, QQ[x_0..x_(#E-1)])
gradf = diff(vars ring f, f)
phi = map(ring f, ring f, gradf);
Phi = rationalMap(phi);
isDominant(phi)
isBirational(phi)
degree Phi

-- 4-cycle plus another node conencted to 3/4 vertices
E = {{0,1}, {1,2}, {2,3}, {0,3}, {0,4}, {1,4}, {2,4}};
G = Graphs$graph E;
edges G
M = matroid G;
f = matroidPolynomial(M, QQ[x_0..x_(#E-1)])
gradf = diff(vars ring f, f)
phi = map(ring f, ring f, gradf);
Phi = rationalMap(phi);
isDominant(phi)
isBirational(phi)
degree Phi