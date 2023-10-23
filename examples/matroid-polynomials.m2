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