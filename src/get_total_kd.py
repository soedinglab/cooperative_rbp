"""
Command-line interface for an easy workflow to calculate the total Kd of RBPs with multiple domains, bases on our analytical calculations and Gillespie simulations.
"""

#Imports
import matplotlib as mpl
import matplotlib.colors as colors

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern} \usepackage{siunitx}']
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Latin Modern Roman']
plt.rcParams['lines.markersize'] = 4


import numpy as np
import argparse


import rbp_model

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('parameter_path', help = "path to a file that stores parameters of the protein")
	parser.add_argument('-s', '--simulate', action='store_true', help = "run gillespie simulations to predict the Kd")
	parser.add_argument('-n', '--replicates', action='store', type = int, default = 1, help = "number of replicate simulations, over which to average, default: 1")
	parser.add_argument('-p', '--plot', action='store_true', help = "plot the result of the simulation in a new window")
	return parser.parse_args()


def load_params(path):
	return rbp_model.get_model_parameters(path)


def simulate_kd(model, num_trajectories, volume):
	results = model.run(False, num_trajectories)

	ordered_names = list(model.sanitized_species_names().keys())

	avg = True
	if avg and num_trajectories > 1:
		trajectories = results[0][:,0], np.average(np.dstack(results)[:,1:], axis = 2), ordered_names
	elif not avg or num_trajectories == 1:
		trajectories = results[0][:,0], results[0][:,1:], ordered_names

	#Gillespy2 (since version 1.3) changes the order of the species in the results array to be sorted by alphabet
	#we need to get the indices of rna, prot and the index range for all bound species
	prot_ind = ordered_names.index('prot')
	rna_ind = ordered_names.index('rna')
	bound_low_ind = ordered_names.index('_' + model.species_names[1]) 
	bound_high_ind = ordered_names.index('_' + model.species_names[-1])

	#print Kd value based on concentrations at the end of the simulation
	print('Total Kd based on the result from the simulation/[M]: ', ((np.mean(rbp_model.pop_to_conc(trajectories[1][-20:-1,prot_ind], volume)) * np.mean(rbp_model.pop_to_conc(trajectories[1][-20:-1,rna_ind], volume))) / (np.mean(rbp_model.pop_to_conc(np.sum(trajectories[1][-20:-1,bound_low_ind:bound_high_ind+1], axis=1), volume)))))

	return trajectories


def analytical_kd(model):
	analytical_kd_result = model.analytical_kd()
	print('Total Kd based on the result from analytical calculations/[M]: ', analytical_kd_result)


if __name__ == '__main__':
	arguments = parse_args()

	params = load_params(arguments.parameter_path)
	model = rbp_model.nxn(*params)

	analytical_kd(model)

	trajectories = None
	if arguments.simulate:
		volume = params[5]
		trajectories = simulate_kd(model, arguments.replicates, volume)

	if arguments.plot:
		if not arguments.simulate or not trajectories:
			raise UserWarning('-p or --plot option given, but no simulation took place. No plot can be shown. Please give this option in conjuction with -s or --simulate option.')
		rbp_model.plot_trajectories(*trajectories)

