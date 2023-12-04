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
adjugate = M -> matrix for i from 1 to numrows(M) list for j from 1 to numcols(M) list ((-1)^(i+j))*(det submatrix'(M, {j-1}, {i-1}))


k = 4 -- very slow for k >= 6
G = Graphs$cycleGraph k
edges G
M = matroid G
f = matroidPolynomial(M, QQ[x_0..x_(k-1)])
eval = map(QQ[x_1..x_(k-1)], ring f, {1, x_1..x_(k-1)})
gradf = diff(vars ring f, f)
gradfprime = eval(gradf)
gradfprime = sub(gradfprime, ring f)
phi = map(ring f, ring f, gradfprime)
Phi = rationalMap(phi)
degree Phi
isDominant(phi)

R = CC[x_0..x_(k-1), y_0..y_(k-1)]
gradf = sub(gradf, ring f)
mat = gradf || (vars R)_{k..(2*k-1)}
mat2 = gradf || matrix {{1,1,1,1,1,1}};
I = minors(2, mat2);
use ring f;
baselocus = ideal flatten entries gradf
baselocus = sub(baselocus, ring f)
I = saturate(I, baselocus);
dim I
J = I + ideal (x_0 - 1);
dim J
sols = solveSystem J_*;
netList sols

T = QQ[t]
zet = (roots (t^4 - 1))_1

hess = diff(transpose(vars ring f), gradf)
coordinates sols_1
eval = map(CC, ring f, {-zet,-zet,3*zet,3*zet,3*zet,3*zet});
H = eval(hess)
det H
img = eval(gradf)
-- x_3 = 1, +++-
--

T = QQ[y_(1,1), y_(1,2), y_(2,2), y_(3,3)]
M = matrix {{y_(1,1), y_(1,2), 0}, {y_(1,2), y_(2,2), -y_(1,2) - y_(2,2)}, {0, -y_(1,2)-y_(2,2), y_(3,3)}}
adjugate M
det M
N = matrix {{x_0 + x_2, -x_2, 0}, {-x_2, x_2 + x_3, -x_3}, {0, -x_3, x_1 + x_3}}
adjugate N