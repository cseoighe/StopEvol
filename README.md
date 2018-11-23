## Selection shapes synonymous stop codon use in mammals


R scripts to implement a model of codon sequence evolution that includes the stop codon, as well as the results of analyses of mammalian coding sequences using this model. 

***

### R Scripts

R scripts to fit and simulate from the model:
* [stopcodon.rscript](https://github.com/cseoighe/StopEvol/blob/master/stopcodon.rscript).  Given a phylogenetic tree and an in-frame sequence alignment (that includes the stop codon as the last position of the alignment), fit the stop-extended codon model.

* [sim_mgf1x4.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_mgf1x4.rscript). Given a phylogenetic tree (with branch lengths, defined as substitutions per codon site) and a textfile containing model parameters, simulate a codon sequence alignment (including stop codon) from the stop-extended codon model, based on the Muse and Gaut codon model and the F1x4 model of codon equilibrium frequencies 

***


### Prerequisites

**R packages:**
* "ape" 
* "expm"

***

### Running the scripts 
**[stopcodon.rscript](https://github.com/cseoighe/StopEvol/blob/master/stopcodon.rscript)**

```
Rscript stopcodon.rscript <treefile.ph> <seqfile.fasta>
```

Where "treefile.ph" contains a phylogenetic tree with branch lengths. Relative branch lengths are treated as fixed and not reoptimized. The script optimizes a scale factor that applies to the overall tree length to allow for differences in branch length units, depending on the evolutionary model.  The input file "seqfile.fasta" contains a codon-aware alignment in FASTA format, which must include the stop codon as the last position of the alignment. Sequences that include a gap or a codon other than a stop codon (i.e. sequences for which the stop codon is not positionally homologous with the last position in the alignment) are excluded.


**[sim_mgf1x4.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_mgf1x4.rscript)**

```
Rscript sim_mgf1x4.rscript treefile parameterfile
```

Where "paramterfile" is a text file containing the following 9 values on a single line: gene_name kappa omega phi pi_A pi_C pi_G pi_T (equilibrium frequencies) alignment_length_in_codons and "treefile.ph" should contain a phylogenetic tree with branch lengths (in units of substitutions per codon site). Using these inputs, a sequence alignment is simulated from a stop-extended codon model, based on the parameterization of Mus & Gaut (1994), with the F1X4 model of codon equilibrium frequencies.


### Output format


**[stopcodon.rscript](https://github.com/cseoighe/StopEvol/blob/master/stopcodon.rscript)**

The script prints optimization progress to screen and, when complete, produces a text file with the same name as <seqfile> but with .stopcodon appended. It contains the following: The maximum log likelihood value obtained and parameter estimates (kappa, omega, treescale, phi), delta_lnL (the difference in maximum log likelihood between a model with phi fixed and with phi a free parameter, and a number (0/1) indicating the successful convergence of the optimization. 


**[sim_mgf1x4.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_mgf1x4.rscript)**

The script produces a sequence alignment file (including the stop codon) in FASTA format, called <gene_name> (from the parameter file) with the extension .sim.mgf1x4.fastas appended.


***


### Example data

**Data to fit the model**

* [ENSG00000111276_CDKN1B.ph](https://github.com/cseoighe/StopEvol/blob/master/ENSG00000111276_CDKN1B.ph)
* [ENSG00000111276_CDKN1B.fasta](https://github.com/cseoighe/StopEvol/blob/master/ENSG00000111276_CDKN1B.fasta)

**Data to simulate from the model**

* [example_parameter_file](https://github.com/cseoighe/StopEvol/blob/master/example_parameter_file)
* [ENSG00000111276_CDKN1B.ph](https://github.com/cseoighe/StopEvol/blob/master/ENSG00000111276_CDKN1B.ph)

***


### Example analysis
```
#Simulate an 100 codon alignment from the model. In the example parameter file phi = 0.1
Rscript sim_mgf1x4.rscript ENSG00000111276_CDKN1B.ph example_parameter_file

#Fit the model
Rscript stopcodon.rscript ENSG00000111276_CDKN1B.ph CDKN1B.sim.mgf1x4.fastas 
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



