loadPackage "RationalMaps"

R = QQ[x,y,z]
phi = map(R, R, {x + y, x - z, z - y})
inverseOfMap(phi)

R = QQ[x,y,z]
psi = map(R, R, {x + y, 2*y + x + z, z + y})
inverseOfMap(psi)