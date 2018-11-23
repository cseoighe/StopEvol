## Fit a mixture model to estimate the proportion of stop codons under purifying selection


R code and data to fit a mixture model to data from the OrthoMaM (v.8) database to estimate the proportion of stop codons under purifying selection pressure

***

### R Scripts

mixture_mo
* [mixture_model.rscript](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/mixture_model.rscript).  Fits the mixture model. 

***

### Prerequisites

**R packages:**
* "ape"
* "expm"

### Input files 

The following files must be located in the same directory as the script:

* [Trees](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Trees). Concatenated trees corresponding to 12,336 orthologue families from the OrthoMaM (v.8) database (only families with at least 20 taxa for which there is a stop codon at the last position of the alignment are included).  
* [Stop_codons](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Stop_codons). The stop codon found at the end of the alignment for each taxon for each orthologue family (the orthologue families appear in the same order as for the concatenated trees). Stop codons are encoded as 1, 2, 3 for UAG, UGA and UAA, respectively (taxa with a gap or a codon other than a stop codon are encoded as NA). 
* [Model_parameters.mgf1x4](https://github.com/cseoighe/StopEvol/blob/master/MixtureModel/Model_parameters.mgf1x4). Maximum likelihood estimates of parameter values, estimated using codonPhyml with the following command: 
codonphyml -i <sequence.phy> -m MG --fmodel F1X4 -t e -f empirical -w DM0 -q -o lr -u <tree.ph>
The file contains the following values: gene_name kappa omega pi_A pi_C pi_G pi_T, where pi_x is the estimated frequency of nucleotide x.

***

### Running the script

Note: This script can take several hours to run for the data provided

```
Rscript mixture_model.rscript 
```


### Authors

* **Cathal Seoighe**
* **Stephen J. Kiniry**
* **Andrew Peters**
* **Pavel V. Baranov**
* **Haixuan Yang**

***

### License

This project is licensed  - see the [LICENSE](LICENSE) file for details

***





