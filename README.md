# Scripts for simulation and plotting of thermodynamic model for cooperative RNA-binding of multidomain proteins

 [![License](https://img.shields.io/github/license/soedinglab/cooperative_rbp.svg)](https://choosealicense.com/licenses/gpl-3.0/)
 [![Issues](https://img.shields.io/github/issues/soedinglab/cooperative_rbp.svg)](https://github.com/soedinglab/cooperative_rbp/issues)

##  Publication
[Stitzinger S.H., Sohrabi-Jahromi S., and Söding J. (2021) Cooperativity boosts affinity and specificity of proteins with multiple RNA-binding domains. NAR Genomics & Bioinformatics](https://doi.org/10.1093/nargab/lqad057).

## Abstract
Numerous cellular processes rely on the binding of proteins with high affinity to specific sets of RNAs. Yet most RNA-binding domains display low specificity and affinity in comparison to DNA-binding domains. The best binding motif is typically only enriched by less than a factor 10 in high-throughput RNA SELEX or RNA bind-n-seq measurements. Here, we provide insight into how cooperative binding of multiple domains in RNA-binding proteins (RBPs) can boost their effective affinity and specificity orders of magnitude higher than their individual domains. We present a thermodynamic model to calculate the effective binding affinity (avidity) for idealized, sequence-specific RBPs with any number of RBDs given the affinities of their isolated domains. For seven proteins in which affinities for individual domains have been measured, the model predictions are in good agreement with measurements. The model also explains how a two-fold difference in binding site density on RNA can increase protein occupancy ten-fold. It is therefore rationalized that local clusters of binding motifs are the physiological binding targets of multi-domain RBPs.

## Simulation script
`rbp_model.py` implements the setup of the simulations and the calculations based on our analytical approach. The parameters of the proteins are loaded from text files in the folder `examples`. 

## Usage with command-line interface
The prediction of the total *K*<sub>d</sub> can be started from the command line, using the script `get_total_kd.py`.

	Usage:
	  python3 get_total_kd.py [-h] [--simulate] [-n REPLICATES] [--plot] parameter_path

	positional arguments:
	  parameter_path	path to a file that stores parameters of the protein

	optional arguments:
	  -h, --help            
	  -s, --simulate	run gillespie simulations to predict the Kd;
	  			without this option, only the calculations are used
	  -n REPLICATES, --replicates REPLICATES
				number of replicate simulations, over which to
				average, default: 1
	  -p, --plot		plot the trajectory of the simulation
	  			and save the plot to PDF, next to the parameter file

For example, to run a simulation of the protein hnRNP A1, averaged over 5 trajectories, and to show a plot of the trajectory after the simulation we can use
	
	python3 get_total_kd.py ../examples/hnrnp_a1.csv -s -n 5 -p

## Plotting script
To plot the figures in our paper the following functions from the script `visualizations.py` were used.

**Figure 1** Many RBPs have more than one domain

	rbp_distribution()
	oligomer_distribution()

**Figure 3C** Effective concentration *c*<sub>12</sub> of domain 2 at RNA site 1, when at least one RNA site is already bound

	ceff_L()
	
**Figure 4** Comparison of experimental *K*<sub>d</sub>s to predictions from our model are in good agreement

	example_overview()

**Figure 5A** Dissociation constants decrease exponentially with the number of binding domains

	kd_tot_N()

**Figure 5B** Individual domains contribute to the total *K*<sub>d</sub> up to a threshold in the individual *K*<sub>d</sub>

	kd_tot_k3()

**Figure 5C** Dissociation constants decrease with the binding site density on the RNA

	kd_motif_density()

**Figure 5D** Relative occupancy as a function of binding site density on the RNA

	occupancy_linker_length()

## Dependencies
The scripts require python3 with the following libraries:
- `numpy`
- `scipy`
- `matplotlib`
- `gillespy2`
