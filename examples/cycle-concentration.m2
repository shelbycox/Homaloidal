-- concentration model for a 4-cycle
-- ML deg > 1, therefore, we expect f is not homaloidal
R = QQ[k_(1,1), k_(1,2), k_(1,4), k_(2,2), k_(2,3), k_(3,3), k_(3,4), k_(4,4)]
K = matrix {{k_(1,1), k_(1,2), 0, k_(1,4)}, {k_(1,2), k_(2,2), k_(2,3), 0}, {0, k_(2,3), k_(3,3), k_(3,4)}, {k_(1,4), 0, k_(3,4), k_(4,4)}}
f = det K
F = factor f;
#F
degree f

v = vars R
V = v ** transpose v;
H = diff(V, f);
h = det H -- runs quickly, non-zero
I = ideal h;
hFactors = factor h;
#hFactors
eval = map(QQ, R, random(QQ^1, QQ^8))
det eval(H)

loadPackage("RationalMaps")
inverseOfMap(map(R, R, diff(v, f))) -- the map is not birational