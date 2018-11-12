## Selection shapes synonymous stop codon use in mammals

Phylogenetic models of the evolution of protein-coding sequences can provide insights into the selection pressures that have shaped them. In the application of these models synonymous nucleotide substitutions, which do not alter the encoded amino acid, are often assumed to have limited functional consequences and used as a proxy for the neutral rate of evolution. The ratio of nonsynonymous to synonymous substitution rates is then used to categorize the selective regime that applies to the protein (e.g. purifying selection, neutral evolution, diversifying selection). Here, we extend these models to explore the extent of purifying selection acting on substitutions between synonymous stop codons. Our results suggest that the preference for UGA stop codons
in many multicellular eukaryotes is selective rather than mutational in origin.

***

### R Scripts

There are two different forms of R script. 

* **First** method, generates maximum likelihood estimates of parameters using a stop-codon-extended model given a phylogenetic tree and a codon-aware alignment ([stopcodon.rscript](https://github.com/cseoighe/StopEvol/blob/master/stopcodon.rscript)).
* **Second** method, uses the stop-codon model to generate an alignment when given a phylogenetic tree and a parameter file ([sim_gyempirical.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_gyempirical.rscript)) and ([sim_mgf1x4.rscript](https://github.com/cseoighe/StopEvol/blob/master/sim_mgf1x4.rscript))

***

### Getting Started

The R script [stopcodon.rscript](https://github.com/cseoighe/StopEvol/blob/master/stopcodon.rscript) requires as input, a codon-aware sequence alignment file in FASTA format and a phylogenetic tree file expaining the relatedness of the alignment.

The tree file should include branch lengths. Note that relative branch lengths are treated as fixed.
The sequence file should contain a codon-aware alignment in fasta format and must include the 
stop codon as the last position of the alignment. Sequences that include a gap or a codon other than
a stop codon (i.e. sequences for which the stop codon is not positionally homologous with the last 
position in the alignment) are excluded.

***

### Prerequisites

* Require tree file (".ph") and a codon-aware sequence alignment (".fasta"). Sample files can be simultaneously obtained from the [OrthoMam](http://www.orthomam.univ-montp2.fr/orthomam/html/) database which is a database of orthologous mammalian markers. Using the “Browse” tab, the “CDS” option at the top of the page was highlighted. Then the “submit” button at the base of the page was pressed. It is possible to limit the tree to specific species based on clicking the icon that accompanies each image, and chromosome by selecting the correct chromosome. Marker’s can be selected by ticking the requisite box in thefar left column. After selection, the accompanying maximun likelihood tree and nucleotide alignment file can be downloaded together by clicking the download icon at the base of the page.  Although it requires a longer computational time, it is advised to use a greater number of taxa in the alignment as the power to detect purifying selection acting on stop codon use is limited by the number of taxa.
* Require "R" packages "ape" and "expm":

These are then called usng the "require" function in the R script:

```
require(ape)
require(expm)
```

***

### Installing

The "ape" and "expm" packages are installed in the "R" environment using the "install.packages()" command.

```
install.packages("ape")
install.packaages("expm")
```

***

### Running the tests

The code is ran in the Linux command line using the following command:

```
Rscript stopcodon.rscript <treefile.ph> <seqfile.fasta>
```

Where "seqfile.fasta" is a codon-aware alignment of the sequence data and the "treefile.ph" is a maximum likelihood estimate of the tree given the codon-aware alignment.

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

### Authors

* **Cathal Seoighe** - *Initial work* - [StopEvol](https://github.com/cseoighe/StopEvol)
* **Stephen J. Kiniry**
* **Pavel V. Baranov**
* **Haixuan Yang**

***

### License

This project is licensed  - see the [LICENSE](LICENSE) file for details

***

### Acknowledgments
