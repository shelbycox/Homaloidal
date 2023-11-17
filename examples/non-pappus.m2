loadPackage("Matroids")
loadPackage("Cremona")

M = specificMatroid("nonpappus")
R = QQ[x_0..x_8];
f = matroidPolynomial(M, R_*);
gradf = diff(vars R, f);
phi = map(R, R, gradf);
Phi = rationalMap phi;
isDominant(Phi)
degree(Phi)