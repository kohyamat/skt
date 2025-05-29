<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Material+Symbols+Outlined:opsz,wght,FILL,GRAD@24,400,0,0&icon_names=open_in_new" >

# SKT

These files are R codes and data files used for calculating competition model parameters of tree species in tropical forest community.


### Paper

Takashi S. Kohyama, Nanako Shigesada, Kokichi Kawasaki, Matthew D. Potts, Zamah S. Nur Hajar, Tetsuo I. Kohyama, and Douglas Sheil (2025) A simple competition model can predict rainforest tree diversity, species abundance and ecosystem functions, Journal of Ecology, 113, 842&#x2013;855. 
<a href="https://doi.org/10.1111/1365-2745.14485">DOI:10.1111/1365-2745.14485</a>



## Contents

The provided dataset of tropical forest plot is a processed subset of the original dataset used in the model parameterisation for tree populations as demonstrated in our paper. Readers interested in using the original Pasoh 50-hectar plot data for purposes other than reviewing our analysis are advised to contact the [Forest Research Institute Malaysia (FRIM)](https://www.frim.gov.my) and the [Center for Tropical Forest Science-Forest Global Earth Observatory (CTFS-Forest GEO)](https://forestgeo.si.edu), Smithsonian Tropical Research Institute. R codes employ 1-ha for subplot size and local species populations ≥ 10 trees per subplot for model parameter estimation (as baseline setting). One easily change subplot size and/or minimum tree density for local population selection.

* **skt.r** - Main R code for estimating species parameters and predicting a stable equilibrium solution in the competition model.

* **functions.r** - R code defining some functions for calculating turnover rates.

* **data/pfr_sub_100by100.csv.gz** - Processed subset of Pasoh 50-ha plot data, Peninsular Malaysia, for ca. 1990 and ca. 2000 censuses.

  * `species` - species ID, integer
  * `subplot` - subplot ID assigned by dividing the entire area into 100m x 100m, integer
  * `dbh1` - stem diameter (cm) in 1990, numeric
  * `dbh2` - stem diameter (cm) in 2000, numeric
  * `t` - inter-census duration (year), numeric

## License

The data files are licensed under the [Creative Commons Attribution 4.0 International License](https://github.com/kohyamat/p-B/blob/main/LICENSE-CC-BY), and the R codes are licensed under the [MIT License](LICENSE).
