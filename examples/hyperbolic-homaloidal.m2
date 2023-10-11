R = QQ[x_1..x_3]
f = sum for i from 1 to 2 do for j from i+1 to 3 list x_i*x_j
f = x_1*x_2 + x_1*x_3 + x_2*x_3

loadPackage("RationalMaps")
phi = map(R,R,diff(vars R, f))
isBirationalMap(phi)

T = QQ[a_1..a_3,b_1..b_3,c_1..c_3]
S = T[x,y,z]
M = matrix {{ x*a_1 + y*b_1 + z*c_1, a_2*x + b_2*y + c_2*z}, {a_2*x + b_2*y + c_2*z, a_3*x + b_3*y + c_3*z}}
D = det M

netList terms D
coeffDetM = ((coefficients D)_1)_0
coeffDetM_0
F = ideal {coeffDetM_0, coeffDetM_2, coeffDetM_5, coeffDetM_1 - coeffDetM_3, coeffDetM_1 - coeffDetM_4}
dim F
degree F
mingens F