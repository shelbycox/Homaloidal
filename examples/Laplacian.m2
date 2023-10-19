loadPackage("Graphs")
-- loadPackage("submatrix")
n = 4 --number of vertices
E = {{0,1},{1,2},{2,3},{3,0},{0,2}} --edge set (start the vertices from 0)
G = graph E 
R = QQ[x_0..x_(length(E)-1)]
--Laplacian matrix
Lap = map(R^n, R^n, (i,j) -> if i==j then sum(positions(E, e -> member(i,e)), ind->x_ind) else if (member({i,j},E) or member({j,i},E)) then -x_(position(E, e-> e=={i,j} or e=={j,i})) else 0)
submatrix'(Lap, {3}, {3})
--The following det is supposed to be our matroid polynomial corresponding to G.
determinant(submatrix'(Lap, {3}, {3}))


