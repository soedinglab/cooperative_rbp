"""
Module for handling input and output and running the simulation
"""


# Imports
from __future__ import division
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import constants

import gillespy
import rbp_model
import pdb_extract as pdbe

#persistence length
lp = 0.9e-9

#length per Base
lpb = 6e-10


def get_model_parameters(model_file, distance=True, pdb_file=None, canonical_file=None, residues=None):
	"""
	Reads the parameters and intial values, which define the model, from a CSV file. See an example file for the structure. The first column contains the label (n, prot0, rna0, on, off, volume, lp, L) and is seperated from the values by a semicolon ';'. The rows of the file should contain the following parameters:
		n: number of bindings sites on protein and RNA
		prot0, rna0: Concentrations [mol L^-1] of unbound Protein and RNA
		On/Off rate constants for uncooperative binding events (Molar) (list)
		Volume
		d: (optional, can be calculated from coordinates in a PDB file) - Euklidian distance between Protein-Domains, Input 1D array with values seperatet by colons (np.array, dimension: n,n)
		L: Distance along RNA chain between the binding sites (list), IN: no. of nucleotides, OUT: dist in metres

	INPUT
		model_file - str - path to the CSV file
		distance - bool - True if a distance matrix is given in CSV, False to calculate distance matrix from PDB coordinates
		pdb_file (optional)- str - path to the PDB file
		canonical_file (optional) - str - path to a file containing the canonical sequence
		residues (optional) - list, int - list of residue numbers in the canocial sequence of the binding sites

	RETURN
		list of model parameters - (n, prot0, rna0, on, off, volume, lp, d, L)
	"""

	params_dict = {}
	params = []

	with open(model_file, 'r') as f:
		for line in f:
			(key, val) = line.split(';')
			params_dict[key] = val
	f.close()


	# add every parameter from the dictionary and check the input for errors (e.g. missing values...)
	try:
		params.append(int(params_dict['n'])) #n
		params.append(float(params_dict['prot0'])) #prot0
		params.append(float(params_dict['rna0'])) #rna0
		params.append([float(i) for i in params_dict['on'].split(',')]) #on
		params.append([float(i) for i in params_dict['off'].split(',')]) #off
		params.append(float(params_dict['volume'])) #volume
		params.append(lp) #lp

		if distance:
			params.append(np.fromstring(params_dict['d'], dtype = float, sep = ',').reshape(2,2)) #d
		else:
			raise SystemExit('Distance calculation not implemented yet. PLease provide a distance matrix.')
			#params.append(pdbe.res_dist(pdb_file, canonical_file, residues)

		params.append([(int(i)*lpb) for i in params_dict['L'].split(',')]) #L

	except (KeyError, ValueError):
		raise SystemExit('Error in Parameter file. Please check for errors and run again.')

	return params



def init_run_model(params, labels=False, num_trajectories=1, avg=True):
	"""
	Creates an instance of the model, runs the simulation and returns the results.
	INPUT
		params - list - parameters to initialize the model (n, prot0, rna0, on, off, volume, lp, L)
		labels - bool - turn labels in trajectories on/off
		num_trajectories - int - number of simulations to be returned
		avg - bool - True to average the results over all trajectories, False to return all tractectories
	OUTPUT
		list - (time, species counts, species names)
	"""
	
	model = rbp_model.nxn(params)
	results = model.run(show_labels=labels, number_of_trajectories=num_trajectories)
	if avg and num_trajectories > 1:
		return (results[0][:,0], np.average(np.dstack(results)[:,1:], axis = 2), model.species_names)
	elif not avg or num_trajectories == 1:
		return (results[0][:,0], results[0][:,1:], model.species_names)



def plot_trajectories((time, species, names)):
	"""
	Creates plots from the simulated data
	INPUT
		time
		species - np.array - results from simulation
		names - list, str - species names
	"""
	
	names.insert(0, 'rna')
	for i in range(species.shape[1]-1):
		plt.plot(time, species[:,i+1], label=names[i])
	
	plt.xlabel('Time [s]')
	plt.ylabel('No. of species')
	plt.legend()
	plt.tight_layout()
	plt.show()



if __name__ == '__main__':
	params = get_model_parameters('zbp1_params.csv')
	volume = params[5]
	
	trajectories = init_run_model(params, num_trajectories = 25)
	plot_trajectories(trajectories)
	print((rbp_model.pop_to_conc(trajectories[1][-1,0], volume) * rbp_model.pop_to_conc(trajectories[1][-1,1], volume)) / rbp_model.pop_to_conc(trajectories[1][-1,-1], volume))
	pass
