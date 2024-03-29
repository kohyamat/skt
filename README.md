# SKT

These files are R codes and data files used for calculating competition model parameters of tree species in tropical forest community.


### Paper

Takashi S. Kohyama, Nanako Shigesada, Kokichi Kawasaki, Matthew D. Potts, Zamah S. Nur Hajar, Tetsuo I. Kohyama, and Douglas Sheil,
"A simple competition model predicts rainforest tree diversity and relative abundance"

## Contents

The provided dataset of tropical forest plot is a processed subset of the original dataset used in the model parameterisation for tree populations as demonstrated in our paper. Readers interested in using the Pasoh 50-hectar plot data for purposes other than reviewing our analysis are advised to contact the [Forest Research Institute Malaysia (FRIM)](https://www.frim.gov.my) and the [Center for Tropical Forest Science-Forest Global Earth Observatory (CTFS-Forest GEO)](https://forestgeo.si.edu), Smithsonian Tropical Research Institute. R codes employ 1-ha for subplot size and local species populations ≥ 10 trees per subplot for model parameter estimation (as baseline setting). One easily change subplot size and/or minimum tree density for local population selection.

* **skt.r** - Main R code for estimating species parameters and predicting a stable equilibrium solution in the competition model.

* **turnover_subplot.r** - R code to calculate turnover rates for each species poulation in a partitioned plot.

* **functions.r** - R code defining some functions for calculating turnover rates.

* **data/pfr.csv** - Pasoh 50-ha plot data, Peninsular Malaysia, for ca. 1990 and ca. 2000 censuses.

  * `Cd` - species code
  * `x, y` - coordinates (m) of stem location
  * `dbh1` - stem diameter (cm) in 1990
  * `dbh2` - stem diameter (cm) in 2000
  * `t` - inter-census duration (year)

* **data/Pasoh_species.csv** - Species code "Cd" in pfr.csv is related to species by "Family", "Genus", and "Species" (scientific name)

* **output/skt_list.csv** - The table of model parameters, with observed and predicted population data (this file will be generated by "skt.r")


## License

The data files are licensed under the [Creative Commons Attribution 4.0 International License](https://github.com/kohyamat/p-B/blob/main/LICENSE-CC-BY), and the R codes are licensed under the [MIT License](LICENSE).
