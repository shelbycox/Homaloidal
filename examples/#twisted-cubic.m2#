R = QQ[X,Y,Z,W][x,y,z,w]
I = ideal {y^2 - x*z, z^2 - y*w, x*w - y*z}
JacI = transpose jacobian I

S = QQ[x,y,z,w,X,Y,Z,W]
JacI = sub(JacI, S)
entries JacI
MM = matrix{{X,Y,Z,W}, (entries JacI)_0, (entries JacI)_1, (entries JacI)_2}
I = sub(I, S)
J = I + minors(3, MM)
J = saturate(J, {x,y,z,w})
eliminate(J, {x,y,z,w})

loadPackage "Resultants"
T = QQ[x,y,z,w]
C = ideal {y^2 - x*z, z^2 - y*w, x*w - y*z}
dualVariety C
