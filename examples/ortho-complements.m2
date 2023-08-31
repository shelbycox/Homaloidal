--==================================================
--- row spans and their orthogonal complements
--==================================================
T = QQ[x,y,z]
M = random(QQ^2, QQ^3)
L = ideal (M*(transpose vars T)) -- linear space is the kernel of M i.e. orthogonal complement of the row span, so the row span is the orthogonal complement of L
rowSpan = ideal ((transpose gens ker M)*(transpose vars T)) -- this is the orthogonal complement of L, the row span of M

dim L == 1 -- orthogonal complement of row span has dim: n - row rank = 1
dim rowSpan == 2 -- row span has dim: row rank = 2
saturate(L + rowSpan, ideal T_*) == ideal 1_T -- the row span and its orthogonal complement intersect only at the origin

P = ideal ((transpose vars T) - (matrix (transpose M)_1)) -- a point in the row span
-- rowSpan contains P
degree localize(rowSpan, P)
dim (rowSpan + P)
degree (rowSpan + P)

Q = ideal ((vars T) - (matrix {{1,1,1}})) -- a point not in the row span
-- rowSpan does not contain Q
degree localize(rowSpan, Q)
dim (rowSpan + Q)
degree (rowSpan + Q)

--=================================================
--- column spans and their orthogonal complement
--=================================================
R = QQ[u,v,w]
N = random(QQ^3, QQ^2)
colSpan = ideal ((transpose gens ker transpose N)*(transpose vars R))
O = ideal ((transpose N)*(transpose vars R))

dim(colSpan) == 2 -- column span had dim: column rank = 2
dim(O) == 1 -- orthogonal complement has dim: kernel rank = 1
saturate(colSpan + O, ideal R_*) == ideal 1_R -- colSpan and its orthogonal complement intersect only at the origin

Q = ideal ((transpose vars R) - (matrix N_1)) -- get a point in the column span
-- colSpan contains Q
degree localize(colSpan, Q)
dim (colSpan + Q)
degree (colSpan + Q)

P = ideal ((vars R) - (matrix {{1,1,1}})) -- a point not in the column span
-- colSpan does not contain P
degree localize(colSpan, P)
dim (colSpan + P)
degree (colSpan + P)