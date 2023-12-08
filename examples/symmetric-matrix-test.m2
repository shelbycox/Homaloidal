n = 3
R = QQ[x_(1,1)..x_(n,n)]
M = matrix for i from 1 to n list for j from 1 to n list x_(min(i,j), max(i,j))
p = diff(vars R, det M)
diff(x_(1,1), det M)
diff(x_(1,2), det M)
diff(x_(2,2), det M)
diff(x_(1,1), det M)
(terms det M)

adjugate = M -> matrix for i from 1 to numrows(M) list for j from 1 to numcols(M) list ((-1)^(i+j))*(det submatrix'(M, {j-1}, {i-1}))

-- This does not come from a chordal graph!!!
R = QQ[x_(1,1), x_(1,2), x_(2,2), x_(3,3)]
M = matrix {{x_(1,1), x_(1,2), 0}, {x_(1,2), x_(2,2), -x_(1,2) - x_(2,2)}, {0, -x_(1,2) - x_(2,2), x_(3,3)}}
A = adjugate M
a = diff(x_(1,1), det M)
b = diff(x_(1,2), det M)
c = diff(x_(2,2), det M)
d = diff(x_(3,3), det M)
e = diff(-x_(1,2) - x_(2,2), det M) -- clearly dependent on b, c!
A_(1,2)
A_(0,1)
A_(1,1)
b - c
f = det M
phi = map(ring f, ring f, diff(vars ring f, f))
vars ring f
isBirational(phi)
diff(vars ring f, f)

-- from 4-cycle plus a chord
R = QQ[y_(1,1), y_(1,2), y_(1,3), y_(2,2), y_(3,3)]
M = matrix {{y_(1,1), y_(1,2), y_(1,3)}, {y_(1,2), y_(2,2), -y_(1,2) -y_(2,2)}, {y_(1,3), -y_(1,2) -y_(2,2), y_(3,3)}}
f = det M
gradf = diff(vars ring f, f);
phi = map(ring f, ring f, gradf);
isBirational(phi)

a = diff(y_(1,2), f) --== 2*A_(0,1) - 2*A_(1,2)
b = diff(y_(2,2), f) --== A_(1,1) - 2*A_(1,2)
c = diff(y_(1,1), f) --== A_(0,0)
d = diff(y_(1,3), f) --== 2*A_(0,2)
e = diff(y_(3,3), f) --== A_(2,2)
I = ideal(a, b, c, d, e)
A_(1,1) % I
A = adjugate M

A_(0,1) + A_(1,2) --- A_(2,2) + A_(0,2)

A_(0,0)*(A_(1,1)*A_(2,2) - A_(1,2)^2) % I
A_(0,1)*(A_(0,1)*A_(2,2) - A_(0,2)*A_(1,2)) % I
A_(1,1)*(A_(0,0)*A_(2,2) - A_(0,2)^2) % I
A_(0,2)*(A_(0,1)*A_(1,2) - A_(0,2)*A_(1,1)) % I

A_(0,1)*A_(1,2) % I

T = QQ[X_(0,0), X_(0,1), X_(1,1), X_(1,2), X_(0,2), X_(2,2)]
J = ideal (X_(0,1) - X_(1,2), X_(1,1) - 2*X_(1,2), X_(0,0), X_(0,2), X_(2,2))
X_(1,2)^2 % J

R = QQ[y_(1,1), y_(1,2), y_(1,3), y_(2,2), y_(2,3), y_(3,3)]/ideal(y_(1,2) + y_(2,2) + y_(2,3))
M = matrix {{y_(1,1), y_(1,2), y_(1,3)}, {y_(1,2), y_(2,2), y_(2,3)}, {y_(1,3), y_(2,3), y_(3,3)}}
f = det M
gradf = diff(vars ring f, f);
adjugate M
vars ring f

R = QQ[y_(1,1), y_(1,2), y_(1,3), y_(2,3), y_(3,3)]
M = matrix {{y_(1,1), y_(1,2), y_(1,3)}, {y_(1,2), -y_(1,2) - y_(2,3), y_(2,3)}, {y_(1,3), y_(2,3), y_(3,3)}}
f = det M
gradf = diff(vars ring f, f)
adjugate M
(gradf_1 - gradf_3)_0
diff(y_(1,2), f)

-- larger chordal graph
R = QQ[y_(1,1), y_(1,3), y_(3,3), y_(3,4), y_(3,5), y_(4,4), y_(4,5), y_(4,6), y_(5,6)]
y55 = -y_(3,5) - y_(4,5) - y_(5,6)
y66 = -y_(4,6) - y_(5,6)
M = matrix {{y_(1,1), y_(1,3), 0, 0, 0},--
            {y_(1,3), y_(3,3), y_(3,4), y_(3,5), 0},--
            {0, y_(3,4), y_(4,4), y_(4,5), y_(4,6)},--
            {0, y_(3,5), y_(4,5), y55, y_(5,6)},--
            {0, 0, y_(4,6), y_(5,6), y66}}
f = det M;
gradf = diff(vars R, f);
A = adjugate M;
phi = map(ring gradf, ring gradf, gradf);
isBirational(phi)

M1 = submatrix(M, {0,1}, {0,1})
f1 = det M1
M2 = submatrix(M, {1,2,3,4}, {1,2,3,4})
f2 = det M2

A_(0,0) == diff(y_(1,1), f)
2*A_(0,1) == diff(y_(1,3), f)
A_(1,1) == diff(y_(3,3), f)
2*A_(1,2) == diff(y_(3,4), f)
A_(2,2) == diff(y_(4,4), f)
(factor A_(0,3))

diff(y_(1,1), f) == det M2
diff(y_(1,1), f1)
A_(0,0) % det M2
(det M % y_(1,1)) % y_(1,3)

T = QQ[a,b,c,d,e]
N = matrix {{a, b, 0}, {b, c, d}, {0, d, e}}
factor det N