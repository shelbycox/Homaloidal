a = 3
b = 2
R = QQ[x_0..x_(a+b+1),t]

dPoly = (d, x) -> sum for i from 0 to d list x_i*t^i
Sabstar = (a, b, R, t) -> resultant(dPoly(a, (R_*)_{0..a}), dPoly(b, (R_*)_{a+1..a+b+1}), t)

f1 = dPoly(a, R_*)
f2 = dPoly(b, (R_*)_{a+1..a+b+1})
g = Sabstar(2, 4, R, t) -- gives us S(a,b)^*
Delta = diff((vars R)_{0..(a+b+1)}, g)

T = QQ[x_0..x_(a+b+1)]
phi = map(T, T, sub(Delta, T));

loadPackage("RationalMaps")
isBirationalMap(phi)

S = QQ[x_0..x_6,y_0..y_6, Degrees=>{{1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}}]
M = Delta || (vars S)_{7..13}