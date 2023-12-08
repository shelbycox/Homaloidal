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


-- Finds the covariance model in a different way.
Model2 = (E,removedvertex) -> (
    SubLap := submatrix'(Lap(E), {removedvertex}, {removedvertex});
    G := Graphs$graph E;
    n := length(vertices G);
    S := {};
    for p in toList((0,0)..(n-2,n-2)) do (if p_0<=p_1 then S = S | {p});
    T := QQ[x_0..x_(length(E)-1), for i in S list s_i];
    SubLap = sub(SubLap, T);
    M := map(T^(n-1), n-1, (i,j) -> if i<=j then s_(i,j) else s_(j,i));
    I := map(T^(n-1), n-1, (i,j) -> if i==j then 1 else 0);
    J := ideal (M * SubLap - I);
    eliminate(toList (x_0..x_(length(E)-1)), J)
)


G = Graphs$graph E;
n = length(vertices G);
S = {};
for p in toList((0,0)..(n-2,n-2)) do (if p_0<=p_1 then S = S | {p})
T = QQ[for i in S list s_i]
M = map(T^4,4,(i,j) -> if i<=j then s_(i,j) else s_(j,i))
f = det submatrix'(M,{1},{0}) - det submatrix'(M,{1},{1}) + det submatrix'(M,{1},{2}) - det submatrix'(M,{1},{3})
g = det submatrix'(M,{3},{0}) - det submatrix'(M,{3},{1}) + det submatrix'(M,{3},{2}) - det submatrix'(M,{3},{3})
h = det submatrix'(M, {3},{0})
isSubset(ideal (f,g,h), J)
saturate(ideal(f,g,h), det M) == J

M1 = map(T^(n-1), n-1, (i,j) -> if i==1 then 1 else if i==n-2 then 1 else if i<=j then s_(i,j) else s_(j,i))
M2 = matrix{{s_(1,1),s_(1,2),s_(1,3)},{s_(1,2),s_(2,2),s_(2,3)},{1,1,1}}
M3 = matrix{{s_(0,1),s_(0,2),s_(0,3)},{s_(1,1),s_(1,2),s_(1,3)},{1,1,1}}
M4 = map(T^(n-1), n-1, (i,j) -> if i==n-3 then 1 else if i==n-2 then 1 else if i<=j then s_(i,j) else s_(j,i))
M5 = map(T^(n-2), n-2, (i,j) -> if i==1 then 1 else if i<=j then s_(i,j) else s_(j,i))
J = sub(Model2(E,3), T)
sub(minors(3,M1), T) + sub(ideal {det M})
isSubset(sub(minors(3,M1), T), J)
isSubset(J,sub(minors(3,M1), T))
sub(minors(3,M1),T)+sub(minors(3,M4), T) == J
sub(minors(3,M1),T)+sub(ideal {det M2},T)+sub(ideal det M3, T) == J
isSubset(sub(minors(3,M4), T), J)
minors(3,M)



E = {{0,1},{1,2},{2,3},{3,0},{0,2},{1,4},{2,4}} --edge set (start the vertices from 0)
-- E = Graphs$edges Graphs$completeGraph 4
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


hes = (n) -> (
 map(QQ^n,QQ^n,(i,j) -> if i==j then 0 else 1 )
)
det hes(6)


E = {{1,2},{2,3},{3,4},{1,4},{0,1},{0,2},{0,3},{0,4}}
Model2(E,0)


