## Selection shapes synonymous stop codon use in mammals


Phylogenetic models of the evolution of protein-coding sequences can provide insights into the selection pressures that have shaped them. In the application of these models synonymous nucleotide substitutions, which do not alter the encoded amino acid, are often assumed to have limited functional consequences and used as a proxy for the neutral rate of evolution. The ratio of nonsynonymous to synonymous substitution rates is then used to categorize the selective regime that applies to the protein (e.g. purifying selection, neutral evolution, diversifying selection). Here, we extend these models to explore the extent of purifying selection acting on substitutions between synonymous stop codons. 

***

### R Scripts

We implemeted the stop-extended codon model in R and provide the following R scripts to fit and simulate from the model:
* Given a phylogenetic tree and an in-frame sequence alignment (that includes the stop codon as the last position of the alignment), obtain maximum likelihood estimates of model parameters
([stopcodon.rscript](https://github.com/cseoighe/StopEvol/blob/master/stopcodon.rscript)).
* Simulate a codon sequence alignment (including stop codon) from the stop-extended codon model, built on the Muse and Gaut codon model and the F1x4 model of codon equilibrium frequencies ([sim_mgf1x4.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_mgf1x4.rscript)). The script requires a phylogenetic tree (with branch lengths, defined as substitutions per codon site) and a textfile containing the parameters of the model,


***

### Getting Started

The R script [stopcodon.rscript](https://github.com/cseoighe/StopEvol/blob/master/stopcodon.rscript) requires a codon-aware sequence alignment file in FASTA format and a phylogenetic tree, with branch lengths, estimated from a cognate standard codon model (the default model implemented in the R script is the Muse and Gaut (1994) model with the F1x4 model of codon equilibrium frequencies).  

Note that relative branch lengths are treated as fixed and not reoptimized. The script optimizes a scale factor that applies to the overall tree length.

The sequence file should contain a codon-aware alignment in FASTA format and must include the 
stop codon as the last position of the alignment. Sequences that include a gap or a codon other than
a stop codon (i.e. sequences for which the stop codon is not positionally homologous with the last 
position in the alignment) are excluded.

The R script [sim_mgf1x4.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_mgf1x4.rscript) requires a paramter file and a phylogenetic tree file as arguments. The parameter file is a textfile containing the following 8 values on a single line: gene_name kappa omega pi_A pi_C pi_G pi_T required_alignment_length_in_codons. The treefile should contain branch lengths (in units of substitutions per codon site). Using these inputs, a sequence alignment is simulated from the stop-extended codon model, based on MG94 X HKY85 X F1x4.
***

### Prerequisites

**All Scripts:**
* Require "R" packages "ape" and "expm":

***

### Information

The following sections are broken down into:

* **Running the tests**
* **Break down into end to end tests**

For the following scripts 
* **[stopcodon.rscript](https://github.com/cseoighe/StopEvol/blob/master/stopcodon.rscript)**
* **[sim_gyempirical.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_gyempirical.rscript)**
* **[sim_mgf1x4.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_mgf1x4.rscript)**

***

### Running the tests
**[stopcodon.rscript](https://github.com/cseoighe/StopEvol/blob/master/stopcodon.rscript)**

The code is run on the command line using the following command:

```
Rscript stopcodon.rscript <treefile.ph> <seqfile.fasta>
```

Where "treefile.ph" contains a phylogenetic tree with branch lengths and "seqfile.fasta" is a codon-aware alignment.

***


### Break down into end to end tests

Ensure all data is in the required directory. Then enter this directory using the "cd" command in the linux command line followed by the path to the data directory.

```
cd /home/YourDataDirectory
```
Run the R script using the Rscript command.

```
Rscript stopcodon.rscript <treefile.ph> <seqfile.fasta>
```

**First Optimisation:** Four numbers will then appear in the terminal window. These iterively optimise the parameters: kappa, omega and treescale factor respectively, with initial starting values of 2 (kappa), 0.2 (omega) and 1 (tresscale). Kappa and omega are rate heterogeneity parameters that model transitions and non-synonymous substitutions respectively. The treescale scaling factor parameter was required to scale the branch lengths of the tree up or down, to account for the use of a different rate matrix in the 64x64 stop-codon extended model when compared with the general 61x61 codon model. The fourth and final value is the maximum likelihood value. These values are calculated letting phi, or the parameter which models mutation rate between alternative stop codons relative to the rate of synonymous substitutions between sense codons, equal to 1.

```
m0.out = optim(c(2,0.2,1),lik_fun,treefile=tree_file,seqfile=seq_file,phifixed=1,control=list(fnscale=-1))
```

**Second Optimisation:** After this optimisation procedure has converged, another optimisation procedure follows, this time allowing the phi variable to be calculated. Five numbers are optimised that appear iteravely in the terminal window:  kappa, omega, treescale, phi with initial values to be optimised: 2, 0.2, 1 and 1 respectively. The fifth and final value to appear is the maximum likelihood value.

```
m1.out = optim(c(2,0.2,1,1),lik_fun,treefile=tree_file,seqfile=seq_file,control=list(fnscale=-1))
```

**Output:** A file called stopcodon.rscript.out, which contains the following: The maximum log likelihood of the second optimisation; ML estimate of kappa of the second optimisation; ML estimate of omega of the second optimisation; ML estimate of the treescaling parameterof the second optimisation; ML estimate of phi delta_lnL (difference in log likelihood compared to a model with phi=1 calculated by getting the difference in maximum likelihood valies between the 2 optimisation procedures; convergence of optimizer (0 = success,1 = failure). 

```
write(c(m1.out$value,m1.out$par,diff,m1.out$convergence),out_file,ncol=10)
```
***

### Running the tests
**[sim_gyempirical.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_gyempirical.rscript)**

The code is run on the command line using the following command:

```
Rscript sim_gyempirical.rscript treefile parameterfile
```

Where "paramterfile" is a text file containing the following 4 values on a single line: gene_name kappa omega alignment_length_in_codons and the "treefile.ph" is a maximum likelihood estimate of the tree given the codon-aware alignment.

***


### Break down into end to end tests

Ensure all data is in the required directory. Then enter this directory using the "cd" command in the linux command line followed by the path to the data directory.

```
cd /home/YourDataDirectory
```
Run the R script using the Rscript command.

```
Rscript sim_gyempirical.rscript treefile parameterfile
```

**Output:** A file called <gene_name>.sim.gyempirical.fastas containing a simulated coding sequence alignment, including
stop codon.

```
<gene_name>.sim.gyempirical.fastas
```

### Running the tests
**[sim_mgf1x4.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_mgf1x4.rscript)**

The code is run on the command line using the following command:

```
Rscript sim_mgf1x4.rscript treefile parameterfile
```

Where "paramterfile" is a text file containing the following 8 values on a single line: gene_name kappa omega pi_A pi_C pi_G pi_T (equilibrium frequencies) alignment_length_in_codons and the "treefile.ph" is a maximum likelihood estimate of the tree given the codon-aware alignment.

***


### Break down into end to end tests

Ensure all data is in the required directory. Then enter this directory using the "cd" command in the linux command line followed by the path to the data directory.

```
cd /home/YourDataDirectory
```
Run the R script using the Rscript command.

```
Rscript sim_mgf1x4.rscript treefile parameterfile
```

**Output:** A file called <gene_name>.sim.mgf1x4.fastas containing a simulated coding sequence alignment, including
stop codon.

```
<gene_name>.sim.mgf1x4.fastas
```
***
### Authors

* **Cathal Seoighe**
* **Stephen J. Kiniry**
* **Pavel V. Baranov**
* **Haixuan Yang**

***

### License

This project is licensed  - see the [LICENSE](LICENSE) file for details

***

### Acknowledgments
* **Andrew Peters (PhD)**

### Test 
![alt text](https://github.com/cseoighe/StopEvol/blob/master/Sim1.png)
