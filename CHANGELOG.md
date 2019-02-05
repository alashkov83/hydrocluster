Ver. 0.2.0
----------
  
* Added rekkergroup property table and support it in all modules
* Removed dbcv scoring function and added s_dbw index
* Added choice for solutions for cluster analysis on GUI based on index values and index extrema
* Added display other solutions for cluster analysis in console
* Added setting - modification interpoint distances, instead clusterization points weights. moddist(u, w) = dist(u, w)/(w u)/2)), where w and u - weighting coefficients of points
* Command line interface of 'testis' has been changed to options on hjson based config file
* Added display information about protein in console and GUI
* Added property tables based on negative or positive charged groups
* Added others noise processing algorithm
* Precomputed sparse neighbors graph was makes problems and was removed. See scikit-learn/scikit-learn#12105
* Corrected pKa on current literature values
* Corrected LIS to LYS in calc_abs_charge function
* Nanodroplet table was normalized to Alanine
* Added hyperanalis and other console utils for scanning solutions
* Added help menu and extended "About" menu for system information
* Updated README.md
* Screenshots was changed
* Publications have been sorted
* Improved GUI
* Minor bugs was fixed and functional was improved


Ver. 0.1.0
----------
First release.
