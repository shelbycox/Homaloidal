loadPackage("Cremona")
R = QQ[y_1..y_5]
M = matrix {{y_1, y_2, y_6}, {y_2, y_3, y_4}, {y_6, y_4, y_5}}
p = det M
gradp = diff(vars R, p)
adjugate = M -> matrix for i from 1 to numrows(M) list for j from 1 to numcols(M) list ((-1)^(i+j))*(det submatrix'(M, {j-1}, {i-1}))
adjM = adjugate(M)
transpose adjM
gradp
phi = map(R, R, gradp)
inverseMap(phi)

N = matrix {{y_1, y_2, y_5}, {y_2, y_3, y_4}, {y_5, y_4, 0}}
adjugate N
q = det N
gradq = diff(vars R, q)

T = QQ[x_1..x_8]
O = matrix {{x_1, 0, x_2}, {x_3, x_4, x_5}, {x_6, x_7, x_8}}
r = det O
gradr = diff(vars T, r)
adjugate 
phir = map(T, T, gradr)
inverseMap(phir)

R = QQ[y_1..y_10]
M = matrix {{y_1, y_2, y_3, y_4}, 
            {y_2, y_5, y_6, y_7}, 
            {y_3, y_6, y_8, y_9},
            {y_4, y_7, y_9, y_10}}
transpose M
p = det M
gradp = diff(vars R, p)
adjugate M
gradp_1