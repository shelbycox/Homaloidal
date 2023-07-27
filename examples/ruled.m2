R = QQ[X,Y,Z,W][x,y,z,w]
I = ideal (x^2 + y^2 - z^2)
JacI = transpose jacobian I

S = QQ[x,y,z,w,X,Y,Z,W]
I = sub(I,S)
JacI = sub(JacI, S)
MM = matrix{{X,Y,Z,W}, (entries JacI)_0}
J = I + minors(2, MM)
J = saturate(J, ideal(x,y,z))
eliminate(J, {x,y,z,w})

loadPackage "Resultants"
T = QQ[x,y,z,w]
C = ideal (x^2 + y^2 - z^2)
dualVariety C
