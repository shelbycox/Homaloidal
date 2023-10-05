R = QQ[x_0..x_6,t]
f1 = x_0 + x_1*t + x_2*t^2 + x_3*t^3
f2 = x_4 + x_5*t + x_6*t^2
g = resultant(f1,f2,t) -- gives us S(a,b)^*
Delta = diff((vars R)_{0..6}, g)

S = QQ[x_0..x_6,y_0..y_6, Degrees=>{{1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {1,0}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,1}}]
M = Delta || (vars S)_{7..13}