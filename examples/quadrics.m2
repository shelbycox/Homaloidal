loadPackage "RationalMaps"

R = QQ[x,y,z]
X = transpose vars R
D = matrix {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}
V = matrix {{36/100, 48/100, -8/10}, {-8/10, 6/10, 0}, {48/100, 64/100, 6/10}}
q = (transpose X)*V*D*(transpose V)*X
J = jacobian q

sub(J_0, {x => 7, z => 1, y => 0})

h = (x + y)*(y + z)*(x + 2*z)
Jh = jacobian h
sub(Jh, {x => -2, y => -1, z => 1})
degree h

texMath D