
--Checks the conjecture: Graph matroid polynomials are homaloidal iff the graph is chordal.
--The conjecture seems to be true for connected graphs.
loadPackage("Matroids")
loadPackage("RationalMaps")
loadPackage("EdgeIdeals")
loadPackage("Cremona")
loadPackage("Graphs")


N = 100
n = 7
for m from 9 to binomial(n,2) do
    (
    print m;
    for j from 1 to N do
        (print j;
        S = QQ[v_1..v_n];
        G = randomGraph(S,m);
        -- print(G);
        EE = EdgeIdeals$edges G;
        E = {};
        for e in EE do 
            (
            E = E | {{index e_0,index e_1}}
            );
        R = QQ[x_0..x_(length(E)-1)];
        -- dismiss "EdgeIdeals";
        -- E = sub(E,R);
        -- loadPackage("Graphs", Reload=>true);
        G = Graphs$graph(E);
        -- loadPackage("Matroids", Reload=>true);
        M = matroid(G);
        B = bases M;
        L = for b in B list product(toList b, j->x_j);
        p = sum(L);
        Jacob = diff(vars R, sub(p,R));
        phi = rationalMap Jacob;
        if (Graphs$isChordal G and not isBirational(phi)) or (not Graphs$isChordal G and isBirational(phi)) then
            print(m,G,E,B)
        )
    )


-- E = {{1,2},{2,3},{3,4},{4,1},{1,3}}

