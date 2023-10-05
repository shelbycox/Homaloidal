loadPackage "DeterminantalRepresentations"
R = RR[x,y,z]
f = x^2 + y^2

g = product gens R
detRep(g, HyperbolicPt => matrix{{1_RR},{1},{1}})
detRep(f)

R = QQ[a_1..a_3,b_1..b_3,c_1..c_3][x,y,z]
M = matrix {{ x*a_1 + y*b_1 + z*c_1, a_2*x + b_2*y + c_2*z}, {a_2*x + b_2*y + c_2*z, a_3*x + b_3*y + c_3*z}}
det M