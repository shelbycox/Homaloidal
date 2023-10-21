loadPackage("Matroids")
loadPackage("Cremona")
loadPackage("Graphs")


-- Input: E, the set of edges of the graph.
-- Outputs the list of coordinates of the inverse map of the Jacobian of the corresponding polynomial.
-- Prints the list of edges. For i=0..|E|-1, x_i corresponds to the i-th edge in the "printed" list (not the input).
inverseOfJacob = (E) -> (
    R := QQ[x_0..x_(length(E)-1)];
    G := Graphs$graph(E);
    M := matroid(G);
    B := bases M;
    L := for b in B list product(toList b, j->x_j);
    p := sum(L);
    Jacob := diff(vars R, sub(p,R));
    phi := rationalMap Jacob;
    M := matrix inverseMap phi;
    print Graphs$edges G;
    for i in 0..numcols(M)-1 list M_(0,i)
)

-- Input: n, the number of vertices of the complete graph.
-- Outputs the list of coordinates of the inverse map of the Jacobian of the corresponding polynomial.
-- Prints the list of edges. For i=0..|E|-1, x_i corresponds to the i-th edge in the "printed" list.
inverseOfJacobComplete = (n) -> (
    E:= Graphs$edges Graphs$completeGraph n;
    inverseOfJacob E
)


-- Example:
E = {{1,2},{2,3},{3,4},{4,1},{1,3}}
inverseOfJacob E

n = 4
inverseOfJacobComplete n