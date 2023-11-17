loadPackage("Graphs")
-- loadPackage("submatrix")

-- E: edge set (start vertices from 0.)
Lap = (E) -> (
    G := Graphs$graph E;
    n := length(vertices G);
    R := QQ[x_0..x_(length(E)-1)];
    --Laplacian matrix
    map(R^n, R^n, (i,j) -> if i==j then sum(positions(E, e -> member(i,e)), ind->x_ind) else if (member({i,j},E) or member({j,i},E)) then -x_(position(E, e-> e=={i,j} or e=={j,i})) else 0)
)

-- E: edge set (start vertices from 0.), removedvertex: vertex removed from Laplacian
-- Gets the covariance model
Model = (E,removedvertex) -> (
    SubLap := submatrix'(Lap(E), {removedvertex}, {removedvertex});
    G := Graphs$graph E;
    n := length(vertices G);
    if removedvertex>=n then (print "error: removed vertex not valid"; break;);
    S := {};
    for p in toList((0,0)..(n-2,n-2)) do (if p_0<=p_1 then S = S | {p});
    T := QQ[x_0..x_(length(E)-1), for i in S list s_i];
    SubLap = sub(SubLap, T);
    M := map(T^2, binomial(n,2), (i,j) -> if i==0 then s_(S_j) else (-1)^(sum toList S_j)*(determinant submatrix'(SubLap, {(S_j)_1}, {(S_j)_0})));
    I := saturate (minors(2,M), ideal({determinant SubLap}));
    J := saturate (I, for i in S list s_i);
    eliminate(toList (x_0..x_(length(E)-1)), J)
)

adjugate = M -> matrix for i from 1 to numrows(M) list for j from 1 to numcols(M) list ((-1)^(i+j))*(det submatrix'(M, {j-1}, {i-1}))



E = {{0,1},{1,2},{0,2},{2,3},{3,4},{2,4}} --edge set (start the vertices from 0)
removedvertex = 3
--The following det is supposed to be our matroid polynomial corresponding to G.
SubLap = submatrix'(Lap(E), {removedvertex}, {removedvertex})
determinant(SubLap)

decompose Model(E,3)

E = {{0,1},{1,2},{2,3},{3,0},{0,2},{1,4},{2,4}}
Lap(E)
removedvertex = 3
SubLap = submatrix'(Lap(E), {removedvertex}, {removedvertex})
det SubLap
inverse SubLap




