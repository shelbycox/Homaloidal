newPackage(
    "Homaloidal",
    Version => "0.1",
    Date => "7/27/23",
    Headline => "Hello World",
    Authors => {{ Name => "Shelby Cox", Email => "spcox@umich.edu", HomePage => "https://websites.umich.edu/~spcox/"}},
    AuxiliaryFiles => false,
    DebuggingMode => false,
    PackageExports => {"Resultants", "Matroids", "Cremona", "Graphs"}
    )

export {"matroidPolynomial"}

-* Code section *-
makeMonomial = (x, v) -> product (for i from 0 to length(entries v)-1 list x_i^(v_i))
matroidPolynomial = (M, x) -> {
    B := basisIndicatorMatrix(M);
    T := for i from 0 to numcols(B)-1 list makeMonomial(x, B_i);
    return sum T
}

-* Documentation section *-
--beginDocumentation()

-* Test section *-
TEST /// -* [insert short title for this test] *-
-- test code and assertions here
-- may have as many TEST sections as needed
///

end--

-* Development section *-
restart
debug needsPackage "Homaloidal"
check "Homaloidal"

uninstallPackage "Homaloidal"
restart
installPackage "Homaloidal"
viewHelp "Homaloidal"
