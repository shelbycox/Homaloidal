-- concentration model for two copies of K_3 with one shared vertex
R = QQ[k_(1,1), k_(1,2), k_(1,3), k_(2,2), k_(2,3), k_(3,3), k_(3,4), k_(3,5), k_(4,4), k_(4,5), k_(5,5)]
K = matrix for i from 1 to 5 list for j from 1 to 5 list if ((i < 3 and j > 3) or (i > 3 and j < 3)) then 0 else k_(min(i,j), max(i,j))
f = det K
F = factor f;
degree f

v = matrix {{k_(1,1), k_(1,2), k_(1,3), k_(2,2), k_(2,3), k_(3,3), k_(3,4), k_(3,5), k_(4,4), k_(4,5), k_(5,5)}}
V = (transpose v) ** v;

H = diff(V, f); -- compute the Hessian
-- evaluate at a point and check if determinant is zero there
eval = map(QQ, R, random(QQ^1, QQ^11))
det eval(H)

h = det H -- takes a long time, but h is non-zero

loadPackage("RationalMaps")
phi = map(R, R, diff(v, f));
inverseOfMap(phi); -- confirming that the jacobian is homaloidal