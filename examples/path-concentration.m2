-- concentration model for two copies of K_3 with one shared vertex
R = QQ[x_0..x_4]
K = matrix {{x_0, x_1, 0}, {x_1, x_2, x_3}, {0, x_3, x_4}}
f = det K
texMath f
F = factor f;
degree f

v = vars R
V = (transpose v) ** v;

H = diff(V, f); -- compute the Hessian
-- evaluate at a point and check if determinant is zero there
eval = map(QQ, R, random(QQ^1, QQ^5))
det eval(H)

h = det H -- takes a long time, but h is non-zero

loadPackage("RationalMaps")
phi = map(R, R, diff(v, f));
psi = inverseOfMap(phi); -- confirming that the jacobian is homaloidal
texMath psi
psi