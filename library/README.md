
# Description

This folder holds the main R functions that all other scripts call. There are four scripts here and two of those are C++ scripts that are called from R.

1. `rcodelib_stableplayer.R` - this script holds several functions that implement the `CREJM variable selection algorithm`.

2. `codelib.cpp` - these are a collection of `C++ functions` that are called by `rcodelib_stableplayer.R`.

3. `rcodelib_stableplayer_postselection.R` - this script holds several functions that implement parameter estimation and prediction under the `CREJM framework`.

4. `codelib_postselection.cpp` - these are `C++ functions` that are called by `rcodelib_stableplayer_postselection.R`.

For any reproducibility analysis, no action is needed from the user as far as these scripts are concerned.
