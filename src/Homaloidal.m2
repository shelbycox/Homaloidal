newPackage(
    "Homaloidal",
    Version => "0.1",
    Date => "7/27/23",
    Headline => "Hello World",
    Authors => {{ Name => "Shelby Cox", Email => "spcox@umich.edu", HomePage => "https://websites.umich.edu/~spcox/"}},
    AuxiliaryFiles => false,
    DebuggingMode => false,
    PackageExports => {"Resultants"}
    )

export {}

-* Code section *-


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
