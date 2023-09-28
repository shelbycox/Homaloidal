S = QQ[A..J]
R = S[x,y,z]
f = A*x^3 + B*y^3 + C*z^3 + D*x^2*y + E*x^2*z + F*y^2*z + G*x*y^2 + H*x*z^2 + I*y*z^2 + J*x*y*z
Jac = diff(matrix {{x,y,z}}, f)
H = diff(transpose(matrix {{x,y,z}}), Jac)
h = det H
coeff = (coefficients h)_1
I = sub(ideal coeff, S) -- variety in P9
dim I -- dim 5 variety in P9
degree I

res I
#(decompose I)