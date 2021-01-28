# Scripts for simulation and plotting of thermodynamic model for cooperative RNA-binding of multidomain proteins

 [![License](https://img.shields.io/github/license/soedinglab/cooperative_rbp.svg)](https://choosealicense.com/licenses/gpl-3.0/)
 [![Issues](https://img.shields.io/github/issues/soedinglab/cooperative_rbp.svg)](https://github.com/soedinglab/cooperative_rbp/issues)

##  Publication
[Stitzinger S.H., Sohrabi-Jahromi S., and Söding J. (2021) Cooperativity boosts affinity and specificity of proteins with multiple RNA-binding domains. bioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.27.428308v1).

## Abstract
Numerous cellular processes rely on the binding of proteins with high affinity to specific sets of RNAs. Yet most RNA binding domains display low specificity and affinity, to the extent that for most RNA-binding domains, the enrichment of the best binding motif measured by HT-SELEX or RNA bind-n-seq is usually below 10-fold, dramatically lower than that of DNA-binding domains. Here, we develop a thermodynamic model to predict the binding affinity for proteins with any number of RNA-binding domains given the affinities of their isolated domains. For the four proteins in which affinities for individual domains have been measured the model predictions are in good agreement with experimental values. The model gives insight into how proteins with multiple RNA-binding domains can reach affinities and specificities orders of magnitude higher than their individual domains. Our results contribute towards resolving the conundrum of missing specificity and affinity of RNA binding proteins and underscore the need for bioinformatic methods that can learn models for multi-domain RNA binding proteins from high-throughput *in-vitro* and *in-vivo* experiments.

## Simulation script
`rbp_model.py` implements the setup of the simulations and the calculations based on our analytical approach. The parameters of the proteins are loaded from text files in the folder `examples`. 


## Plotting script
To plot the figures in our paper the following functions from the script `visualizations.py` were used.

**Figure 1** Many RBPs have more than one domain

	rbp_distribution()
	
**Figure 3** Comparison of experimental *K*<sub>d</sub>s to predictions from our model are in good agreement

	example_overview()

**Figure 4A** Dissociation constants decrease exponentially with the number of binding domains

	kd_tot_N()

**Figure 4B** Individual domains contribute to the total *K*<sub>d</sub> up to a threshold in the individual *K*<sub>d</sub>

	kd_tot_k3()

**Figure 4C** Dissociation constants decrease with the binding site density on the RNA

	kd_motif_density()

**Figure 4D** Relative occupancy as a function of binding site density on the RNA

	occupancy_linker_length()


## Dependencies
The scripts require python3 with the following libraries:
- `numpy`
- `scipy`
- `matplotlib`
- `gillespy2`
