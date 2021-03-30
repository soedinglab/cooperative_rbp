# Scripts for simulation and plotting of thermodynamic model for cooperative RNA-binding of multidomain proteins

##  Publication
[Stitzinger S.H., Sohrabi-Jahromi S., and Söding J. (2021) Cooperativity boosts affinity and specificity of proteins with multiple RNA-binding domains. bioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.27.428308v1).

## Abstract
Numerous cellular processes rely on the binding of proteins with high affinity to specific sets of RNAs. Yet most RNA binding domains display low specificity and affinity, to the extent that for most RNA-binding domains, the enrichment of the best binding motif measured by HT-SELEX or RNA bind-n-seq is usually below 10-fold, dramatically lower than that of DNA-binding domains. Here, we develop a thermodynamic model to predict the binding affinity for proteins with any number of RNA-binding domains given the affinities of their isolated domains. For the four proteins in which affinities for individual domains have been measured the model predictions are in good agreement with experimental values. The model gives insight into how proteins with multiple RNA-binding domains can reach affinities and specificities orders of magnitude higher than their individual domains. Our results contribute towards resolving the conundrum of missing specificity and affinity of RNA binding proteins and underscore the need for bioinformatic methods that can learn models for multi-domain RNA binding proteins from high-throughput *in-vitro* and *in-vivo* experiments.

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
				average, default: 10
	  -p, --plot		plot the trajectory of the simulation in a new window

For example, to run a simulation of the protein hnRNP A1, averaged over 5 trajectories, and to show a plot of the trajectory after the simulation we can use
	
	python3 get_total_kd.py ../example/hnrnp_a1.csv -s -n 5 -p

## Plotting script
To plot the figures in our paper the following functions from the script `visualizations.py` were used.

**Figure 1** Many RBPs have more than one domain

	rbp_distribution()
	oligomer_distribution()

**Figure 3A** Dissociation constants decrease exponentially with the number of binding domains

	kd_tot_N()

**Figure 3B** Individual domains contribute to the total *K*<sub>d</sub> up to a threshold in the individual *K*<sub>d</sub>

	kd_tot_k3()

**Figure 3C** Dissociation constants decrease with the binding site density on the RNA

	kd_motif_density()

**Figure 3D** Relative occupancy as a function of binding site density on the RNA

	occupancy_linker_length()

**Figure 4** Comparison of experimental *K*<sub>d</sub>s to predictions from our model are in good agreement

	example_overview()

## Dependencies
The scripts require python3 with the following libraries:
- `numpy`
- `scipy`
- `matplotlib`
- `gillespy2`
