## Perform a bootstrap over alignments to assess uncertainty in the estimates of p and phi in the mixture model.


R code and data (same data as the mixture model uses) in order to perform a bootstrap over alignments to assess uncertainty in the estimates of p and phi in the mixture model.

***

### R Scripts

[bootstrap.rscript](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Bootstraps/bootstrap.rscript)

***

### Prerequisites

**R packages:**
* "ape"
* "expm"

### Input files 

The following files (used in the mixture model) must be located in the parent directory of the script:

* [Trees](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Trees). Concatenated trees corresponding to 12,336 orthologue families from the OrthoMaM (v.8) database (only families with at least 20 taxa for which there is a stop codon at the last position of the alignment are included).  
* [Stop_codons](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Stop_codons). The stop codon found at the end of the alignment for each taxon for each orthologue family (the orthologue families appear in the same order as for the concatenated trees). Stop codons are encoded as 1, 2, 3 for UAG, UGA and UAA, respectively (taxa with a gap or a codon other than a stop codon are encoded as NA). 
* [Model_parameters.mgf1x4](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Model_parameters.mgf1x4). Maximum likelihood estimates of parameter values, estimated using codonPhyml with the following command: 
codonphyml -i <sequence.phy> -m MG --fmodel F1X4 -t e -f empirical -w DM0 -q -o lr -u <tree.ph>
The file contains the following values: gene_name kappa omega pi_A pi_C pi_G pi_T, where pi_x is the estimated frequency of nucleotide x.

***

### Running the script

Note: This script can take a day to run for the data provided

```
Rscript bootstrap.rscript
```


***

### Output 

Produces 4 RData files - "mixture_model_se_bs1000[all|UGA|UAG|UAA].RData":

* [All](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Bootstraps/mixture_model_se_bs1000all.RData)
* [UGA](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Bootstraps/mixture_model_se_bs1000UGA.RData)
* [UAG](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Bootstraps/mixture_model_se_bs1000UAG.RData)
* [UAA](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Bootstraps/mixture_model_se_bs1000UAA.RData)

Each contains 1000 bootstrap estimates for phis and weights. 





