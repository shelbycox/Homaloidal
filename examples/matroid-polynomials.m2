loadPackage("Matroids")
loadPackage("RationalMaps")
loadPackage("Graphs")
M = matroid({a,b,c,d,e}, {{a,b,d}, {c,d,e}}, EntryMode=>"nonbases")
M = matroid({a,b,c,d},{}, EntryMode => "nonbases")
U23 = uniformMatroid(2, 3)
U24 = uniformMatroid(2, 4)
B = basisIndicatorMatrix U24

makeMonomial = (x, v) -> product (for i from 0 to length(entries v)-1 list x_i^(v_i))
matroidPolynomial = (M, x) -> {
    B := basisIndicatorMatrix(M);
    T := for i from 0 to numcols(B)-1 list makeMonomial(x, B_i);
    return sum T
}

graph ({{1,2},{2,3},{3,4}})

M = specificMatroid("vamos")
R = QQ[x_0..x_7]
f = matroidPolynomial(M, x)
#(factor f)
phi = map(R, R, diff(vars R, f))
isBirationalMap(phi)

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