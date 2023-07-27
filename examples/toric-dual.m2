R = QQ[X,Y,Z,W][x,y,z,w]
f = x*w - y*z
M = transpose jacobian f

S = QQ[x,y,z,w,X,Y,Z,W]
M = sub(M,S) 
f = sub(f,S)
MM = matrix{{X,Y,Z,W}, (entries M)_0}
J = ideal f + minors(2, MM)
J = saturate(J, {x,y,z,w})
eliminate(J, {x,y,z,w})

loadPackage("Resultants", Reload => true)
T = QQ[x,y,z,w]
I = ideal(x*w - y*z)
dualVariety I
