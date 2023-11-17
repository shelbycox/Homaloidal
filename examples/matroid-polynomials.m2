loadPackage("Matroids")
loadPackage("Cremona")
loadPackage("Graphs")
M = matroid({a,b,c,d,e,f}, {{a,b,c,d,e}, {a,b,c,d,f}, {a,b,c,e,f}}, EntryMode => "bases")
M = matroid({a,b,c,d},{{abcd}}, EntryMode => "nonbases")
U23 = uniformMatroid(2, 3)
U34 = uniformMatroid(3, 4)
B = basisIndicatorMatrix U24

makeMonomial = (x, v) -> product (for i from 0 to length(entries v)-1 list x_i^(v_i))
matroidPolynomial = (M, x) -> {
    B := basisIndicatorMatrix(M);
    T := for i from 0 to numcols(B)-1 list makeMonomial(x, B_i);
    return sum T
}

G = Graphs$graph {{1,2},{2,3},{3,4},{1,4},{1,3},{2,4}}
M = matroid G
#(bases M)
edges G

R = QQ[x_0..x_5]
f = matroidPolynomial(M, x)
#(factor f)
phi = map(R, R, diff(vars R, f))
psi = inverseMap(phi)
g = ((matrix psi)_0)_0
factor phi(g)
peek psi
PHI = inverseMap(psi)
PHI*psi
compose(phi, psi)
compose(psi, phi)


--           12 13 14 15 23 24 25 34 35 45
D = matrix {{0, 0, 0, 0, 1, 0, 1, 1, 1, 1}, 
            {0, 1, 0, 1, 0, 0, 0, 1, 1, 1},
            {1, 0, 1, 1, 0, 1, 1, 0, 0, 0},
            {0, 1, 0, 1, 1, 0, 1, 0, 0, 0},
            {1, 1, 1, 0, 1, 1, 0, 0, 0, 0}}
D = sub(D, QQ)
reducedRowEchelonForm(D)

-- YES homaloidal: U2*, 
-- NOT homaloidal: C5, fano, P6
-- ??? homaloidal: vamos

K6 = Graphs$completeGraph(6)
M6 = matroid K6
R6 = QQ[x_0..x_14]
f6 = matroidPolynomial(M6, R6_*);
phi6 = map(R6, R6, diff(vars R6, f6));
psi6 = inverseMap(phi6)

Mp = specificMatroid("nonpappus")
circuits Mp
flats Mp
-- non-pappus matroid is chordal!
Rp = QQ[x_0..x_8]
fp = matroidPolynomial(Mp, Rp_*);
phip = map(Rp, Rp, diff(vars Rp, fp));
isBirational(phip) -- not birational
-- what is the image? is a section birational?
-- regular + chordal matroid is graphic?

-- two triangles glued at a vertex
E = {{0,1}, {1,2}, {0,2}, {2,3}, {3,4}, {2,4}};
G = Graphs$graph E;
M = matroid G;
R = QQ[x_0..x_5];
f = matroidPolynomial(M, R_*);
phi = map(R, R, diff(vars R, f));
factor phi(f)
psi = inverseMap(phi)
g = ((matrix psi)_2)_0
factor g
factor phi(g)
h2 = x_0 + x_1 - x_2
factor phi(h)

-- triangle and an edge
E = {{0,1}, {1,2}, {0,2}, {2,3}}
G = Graphs$graph E
M = matroid G
R = QQ[x_0..x_3]
f = matroidPolynomial(M, R_*)
phi = map(R, R, diff(vars R, f))
psi = inverseMap(phi)
g = ((matrix psi)_0)_0
factor g
factor phi(g)

-- four cycle with a chord
E = {{0,1}, {1,2}, {2,3}, {0,3}, {0,2}}
G = Graphs$graph E
M = matroid G
R = QQ[x_0..x_4]
f = matroidPolynomial(M, R_*)
phi = map(R, R, diff(vars R, f))
Phi = rationalMap(phi)
degree(Phi)
psi = inverseMap(phi)
g = ((matrix psi)_0)_0
factor g
factor phi(g)

-- four cycle
E = {{0,1}, {1,2}, {2,3}, {0,3}}
G = Graphs$graph E
edges G
M = matroid G
R = QQ[x_0..x_3]
f = matroidPolynomial(M, R_*)
phi = map(R, R, diff(vars R, f))
Phi = rationalMap(phi)
degree(Phi)

-- five cycle
E = {{0,1}, {1,2}, {2,3}, {3,4}, {0,4}}
G = Graphs$graph E
edges G
M = matroid G
R = QQ[x_0..x_4]
f = matroidPolynomial(M, R_*)
phi = map(R, R, diff(vars R, f))
Phi = rationalMap(phi)
degree(Phi)

-- complete graph on four vertices
E = {{0,1}, {1,2}, {2,3}, {0,3}, {0,2}, {1,3}}
G = Graphs$graph E
M = matroid G
R = QQ[x_0..x_5]
f = matroidPolynomial(M, R_*)
phi = map(R, R, diff(vars R, f))
psi = inverseMap(phi)
g = ((matrix psi)_0)_0
factor g
factor phi(g)

-- triangle
E = {{0,1}, {1,2}, {0,2}}
G = Graphs$graph E
M = matroid G
R = QQ[x_0..x_2]
f = matroidPolynomial(M, R_*)
phi = map(R, R, diff(vars R, f))
psi = inverseMap(phi)
g = ((matrix psi)_0)_0
factor phi(g)

-- non-pappus
M = specificMatroid("nonpappus")
R = QQ[x_0..x_8];
f = matroidPolynomial(M, R_*);
gradf = diff(vars R, f);
phi = map(R, R, gradf);
Phi = rationalMap phi;
isDominant(Phi)
degree(Phi)