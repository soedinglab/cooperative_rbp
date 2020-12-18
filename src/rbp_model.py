"""
Gillespie model for cooperative binding of a Protein to ssRNA. N binding sites on protein and RNA
Includes a class for Gillespie simulations and a method for calculating the Kd based on the analytical solution.
"""

# Imports
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import constants
import sys

import gillespy2 as gillespy

#review of RNA chain flexibility in Bao, 2016
#these values are the mean of measurements from five studies
#persistence length RNA
lp = 2.7e-9
lp_err = 0.6e-9

#length per base
lpb = 5.5e-10
lpb_err = 0.9e-10

# protein chain flexibility parameters in Zhou, 2001
# persistence length protein
lp_p = 3.04e-10

# length per amino acid
lpaa = 3.8e-10

class nxn(gillespy.Model):
	"""
	Class for defining model parameters, species and reactions. Inherits from the model class of the gillespy library.
	Methods (for external use):
		analytical_kd - calculate the Kd based on the analytical solution
	"""

	def __init__(self, n, prot0, rna0, on, off, volume, lp, d, L, L_p, lp_p, time):
		"""
		Initialize the model. Initial values are passed when creating a class instance:
			n: number of bindings sites on protein and RNA
			prot0, rna0: Concentrations [mol L^-1] of unbound Protein and RNA
			On/Off rate constants for uncooperative binding events (Molar) (list)
			Volume
			lp: Persistence length RNA
			d: Euklidian distance between Protein-Domains (np.array, dimension: n,n)
			L: Distance along RNA chain between the binding sites (list)
			L_p: Distance along flexible linker between binding domains
			lp_p: Persistence length protein
		"""

		# Initial values passed into the class
		self.n = n
		self.prot0 = prot0
		self.rna0 = rna0
		self.on = on
		self.off = off
		self.my_volume = volume
		self.lp = lp
		self.d = d
		self.L = L
		self.L_p = L_p
		self.lp_p = lp_p
		self.lp_err = lp_err
		self.lpb_err = lpb_err


		gillespy.Model.__init__(self, name="nxn")


		# Parameters:
		#-------------
		# convert Molar On rate constants to stochastic rate constants

		self.on_stoch_bi = []

		self.on_stoch_uni = []
		self.off_stoch = []
		for i in range(self.n):
			self.on_stoch_bi.append(gillespy.Parameter(
				name="on_stoch_bi" + str(i),
				expression=(self.on[i])/(self.my_volume*constants.N_A)))
			self.on_stoch_uni.append(gillespy.Parameter(
				name="on_stoch_uni" + str(i),
				expression=self.on[i]))
			self.off_stoch.append(gillespy.Parameter(
				name="off_stoch" + str(i),
				expression=self.off[i]))


		self.add_parameter(self.on_stoch_uni)
		self.add_parameter(self.on_stoch_bi)
		self.add_parameter(self.off_stoch)


		# Species
		#---------

		self.prot = gillespy.Species(name='prot', initial_value=conc_to_pop(self.prot0, self.my_volume))
		self.rna = gillespy.Species(name='rna', initial_value=conc_to_pop(self.rna0, self.my_volume))
		self.add_species([self.prot, self.rna])

		self.species_names = []
		self.species = []
		for i in range(2**(self.n)):
			self.species_names.append(bin(i)[2:].zfill(self.n))
			self.species.append(gillespy.Species(name="_" + self.species_names[-1], initial_value = 0))
		
		self.add_species(self.species)
		self.delete_species('_' + self.species_names[0])


		# Reactions
		#-----------
		self.add_reaction(self.create_reactions())


		self.timespan(np.linspace(0,time[0],time[1]))


	def run(self, labels, num_trajectories):
		"""Overrides the run method from the parent class. May print a warning and then runs the original method.
		"""
		#Print warning for proteins with a flexible linker between binding domains and more than two binding sites. This can not be done easily in the simulation. The analytical result can be used instead.
		if self.simulation_warning == True:
			print('Please note, that simulations are not implemented for proteins containing a flexible linker between binding domains and more than two binding sites. Please use the function analytical_kd to estimate the total kd based on our analytical result.')
		return super().run(show_labels=labels, number_of_trajectories=num_trajectories)



	def gauss_chain(self, x, mu, sig_sq):
		return ((1/(2*constants.pi*sig_sq))**(3/2)*np.exp(-((x-mu)**2)/(2*sig_sq))) #radial distribution function of a worm like chain in the gaussian limit

	def flex_peptide_chain(self, R1, R2, L, L_p):
		sig_sq = (2/3)*self.lp*L + (2/3)*self.lp_p*L_p
		return((1/((2*constants.pi*sig_sq)**(3/2))) * ((np.exp(-(1/(2*sig_sq))*(R1-R2)**2) - np.exp(-(1/(2*sig_sq))*(R1+R2)**2)) / (2*R1*R2/sig_sq)))

	def get_concentration(self, domain1, domain2):
		"""Return the concentration of domain2 at binding site 2 if domain1 is already bound at site 1 and for symmetrie reasons also vice versa.
		INPUT
			domain1,2 - int - index of the binding domains, zero based index
		"""
		L_tot = 0
		d = 0
		L_p_tot = 0
		R1 = 0
		R2 = 0
		c_eff = 0

		#Ensure domain1<domain2 and exchange the order if not, the c_eff equations are symmetric so we are allowed to do this
		if domain1 > domain2:
			tmp = domain1
			domain1 = domain2
			domain2 = tmp

		L_tot = sum(self.L[domain1*2:(domain2*2-1)]) #RNA linker length
		L_p_tot = sum(self.L_p[domain1:domain2]) + np.trace(self.d[domain1+1:domain2, domain1+1:domain2]) #Protein linker length

		if L_p_tot == 0:
			d = self.d[domain1, domain2]
			sig_sq = (2/3) * self.lp * L_tot
			c_eff = self.gauss_chain(d, 0, sig_sq) * 10**(-3) / constants.N_A
		elif L_p_tot !=0:
			R1 = self.d[domain1, domain2]
			R2 = self.d[domain2, domain1]
			c_eff = self.flex_peptide_chain(R1, R2, L_tot, L_p_tot)* 10**(-3) / constants.N_A

		return c_eff


	def create_reactions(self):
		"""
		Permute all possible reactions and determine the type of reaction, return a list of reactions as gillespy reaction-objects
		"""

		reactions = []
		self.simulation_warning = False

		#loop all species
		for s in self.species_names: 
			
			#find empty binding sites per species
			pos_reac = [] #list of all empty binding sites in current s
			for ind, elem in enumerate(s):
				if elem == '0':
					pos_reac.append(ind)

			#loop all empty binding sites per species and determine the possible reaction plus reverse reaction
			for pos_reac_ind, reac_id in enumerate(pos_reac):

				# Print warning at the end of initialization for proteins with a flexible linker between binding domains and more than two binding sites. This can not be done easily in the simulation. The analytical result can be used instead.
				if np.sum(self.L_p) != 0 and self.n > 2:
					self.simulation_warning = True


				#determine the product of reaction, first as string, then the species object
				product_name = list(s)
				product_name[reac_id] = '1'
				product_name = "".join(product_name)
				product_object = self.species[int(product_name, 2)]
				reactant_object = self.species[int(s, 2)]

				#first binding, n free binding sites
				if len(pos_reac) == (self.n):
					#forward reaction
					reactions.append(gillespy.Reaction(
						name='r' + str(len(reactions)+1),
						reactants={self.prot:1, self.rna:1},
						products={product_object:1},
						rate=self.on_stoch_bi[reac_id]))

					#reverse reaction
					reactions.append(gillespy.Reaction(
						name='r' + str(len(reactions)+1),
						reactants={product_object:1},
						products={self.prot:1, self.rna:1},
						rate=self.off_stoch[reac_id]))

				#second and nth binding, n-1 free binding sites
				elif len(pos_reac) < (self.n) and len(pos_reac) != 0:

					#first check if free binding site at reac_id is in between two already occupied binding sites -> is_middle
					is_middle = False
					l = reac_id
					r = reac_id
					r_neighbour = False
					l_neighbour = False
					while (l > 0 or r < (len(s)-1)):
						l -= 1
						r += 1
						if l >= 0:
							if l_neighbour == False and s[l] == '1':
								l_neighbour = True
						if r < (len(s)):
							if r_neighbour == False and s[r] == '1':
								r_neighbour = True

					if r_neighbour == True and l_neighbour == True:
						is_middle = True


					#binding site is not in the middle of two already bound sites
					if not is_middle:
						#find the bound binding site closest to the reac_id'th binding site
						bound_neighbour = np.nan
						l = reac_id
						r = reac_id
						while (l > 0 or r < (len(s)-1)) and np.isnan(bound_neighbour):
							l -= 1
							r += 1
							if l >= 0 and int(s[l]) == 1:
								bound_neighbour = l
							elif r <= (len(s)-1) and int(s[r]) == 1:
								bound_neighbour = r


						#add gaussian chain distr as parameter
						c_eff = self.get_concentration(reac_id, bound_neighbour)
						self.add_parameter(gillespy.Parameter(
							name="chain_distr"+ str(len(reactions)+1),
							expression=(c_eff)))

						#forward reaction, propensity function (example for 01 -> 11): [01]*on_stoch_1_uni*chain_distr
						reactions.append(gillespy.Reaction(
							name='r' + str(len(reactions)+1),
							reactants={reactant_object:1},
							products={product_object:1},
							propensity_function= '_' + str(s) + '*' + 'chain_distr' + str(len(reactions)+1) + '*' + 'on_stoch_uni' + str(reac_id)))

						#reverse reaction
						reactions.append(gillespy.Reaction(
							name='r' + str(len(reactions)+1),
							reactants={product_object:1},
							products={reactant_object:1},
							rate=self.off_stoch[reac_id]))

					#binding site is in the middle of two already bound sites
					if is_middle:
						#find left and right bound neighbours
						l_bound_neighbour = np.nan
						r_bound_neighbour = np.nan
						l = reac_id
						r = reac_id
						while (l > 0 or r < (len(s)-1)) and (np.isnan(l_bound_neighbour) or np.isnan(r_bound_neighbour)):
							l -= 1
							r += 1
							if l >= 0 and s[l] == '1':
								l_bound_neighbour = l
							if r <= (len(s)-1) and s[r] == '1':
								r_bound_neighbour = r

						

						#total distance between binding sites along RNA chain (L_tot) and between binding domains of the protein (d_tot), left and right of the bound site
						l_d_tot = np.nan
						r_d_tot = np.nan
						l_L_tot = np.nan
						r_L_tot = np.nan

						l_d_tot = self.d[l_bound_neighbour, reac_id]
						l_L_tot = sum(self.L[l_bound_neighbour*2:reac_id*2])

						r_d_tot = self.d[r_bound_neighbour, reac_id]
						r_L_tot = sum(self.L[reac_id*2-1:r_bound_neighbour*2])

						#add gaussian chain distr as parameter
						sig_sq_l = (2/3) * self.lp * l_L_tot
						sig_sq_r = (2/3) * self.lp * r_L_tot
						mu = ((l_d_tot + r_d_tot) * sig_sq_l)/(sig_sq_l + sig_sq_r)
						sig_sq = (sig_sq_l * sig_sq_r)/(sig_sq_l + sig_sq_r)
						self.add_parameter(gillespy.Parameter(
							name="chain_distr"+ str(len(reactions)+1),
							expression=(((self.gauss_chain(l_d_tot, mu, sig_sq))*10**(-3))/constants.N_A)))

						#forward reaction, propensity function (example for 101 -> 111): [01]*on_stoch_1*chain_distr
						reactions.append(gillespy.Reaction(
							name='r' + str(len(reactions)+1),
							reactants={reactant_object:1},
							products={product_object:1},
							propensity_function= '_' + str(s) + '*' + 'chain_distr' + str(len(reactions)+1) + '*' + 'on_stoch_uni' + str(reac_id)))

						#reverse reaction
						reactions.append(gillespy.Reaction(
							name='r' + str(len(reactions)+1),
							reactants={product_object:1},
							products={reactant_object:1},
							rate=self.off_stoch[reac_id]))


		return reactions

	
	def analytical_kd(self):
		"""
		Return the Kd calculated from the analytical solution. Class needs to be initialized with model parameters.
		"""
		kd = 0
		# loop all species (binding configurations) and determine kd for theoretical one step reaction from the unbound state to the i'th species (binding configuration), total kd is the  inverse of the sum of all these individual kds
		for s in self.species_names:
			bound_sites = []
			L_tot = 0
			d = 0
			L_p_tot = 0
			R1 = 0
			R2 = 0
			c_eff = []
			prev_ind = np.nan
			for ind, elem in enumerate(s):
				if elem == '1':
					bound_sites.append(ind)
					if not np.isnan(prev_ind): # find pairs of bound sites, then we to calculate c_eff
						c_eff.append(self.get_concentration(prev_ind, ind))

					prev_ind = ind

			kd += np.prod(np.array(self.on)[bound_sites]) * np.prod(np.array(c_eff))

		return (kd**(-1))



def get_model_parameters(model_file):
	"""
	Reads the parameters and intial values, which define the model, from a CSV file. See an example file for the structure. The first column contains the label (n, prot0, rna0, on, off, volume, lp, L, time) and is seperated from the values by a semicolon ';'. The rows of the file should contain the following parameters:
		n: number of bindings sites on protein and RNA
		prot0, rna0: Concentrations [mol L^-1] of unbound Protein and RNA
		On/Off rate constants for uncooperative binding events (Molar) (list) (optional, Kd values for individual domains can be used instead)
		Kd (optional, on/off rate constants can be used instead) - list - Kd values for individual domains
		Volume [L]
		d: Euklidian distance between Protein-Domains, Input 1D array with values seperatet by colons (np.array, dimension: n,n)
		L: Distance along RNA chain between the binding sites (list), IN: no. of nucleotides, OUT: dist in metres
		time (2 values): end point of simulation, number of time steps to store

	INPUT
		model_file - str - path to the CSV file

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
		if 'on' in params_dict and 'off' in params_dict:
			params.append([float(i) for i in params_dict['on'].split(',')]) #on
			params.append([float(i) for i in params_dict['off'].split(',')]) #off
		elif 'Kd' in params_dict:
			params.append([(1/float(i)) for i in params_dict['Kd'].split(',')]) #on
			params.append([1 for i in params_dict['Kd'].split(',')]) #off
		params.append(float(params_dict['volume'])) #volume
		params.append(lp) #persistence length RNA
		params.append(np.fromstring(params_dict['d'], dtype = float, sep = ',').reshape(params[0],params[0])) #d
		params.append([(int(i)*lpb) for i in params_dict['L'].split(',')]) #length of RNA linkers
		params.append([(int(i)*lpaa) for i in params_dict['L_p'].split(',')]) #length of flexible protein linkers, lpaa (length per amino acid)
		params.append(lp_p) #persistence length protein
		params.append([int(i) for i in params_dict['time'].split(',')]) #end time, timesteps

	except (KeyError, ValueError):
		raise SystemExit('Error in Parameter file. Please check for errors and run again.')

	return params



def init_run_print_model(parameter_file, labels=False, num_trajectories=1, avg=True, simulate=True, plot = True):
	"""
	Creates an instance of the model, runs the simulation and returns the results.
	INPUT
		params - list - parameters to initialize the model (n, prot0, rna0, on, off, volume, lp, L)
		labels - bool - turn labels in trajectories on/off
		num_trajectories - int - number of simulations to be returned
		avg - bool - True to average the results over all trajectories, False to return all trajectories
	OUTPUT
		list - (time, species counts, species names)
	"""
	
	params = get_model_parameters(parameter_file)
	volume = params[5]
	model = nxn(*params)
	analytical_kd_result = model.analytical_kd()
	print('Total Kd based on the result from analytical calculations: ', analytical_kd_result)
	results = model.run(labels, num_trajectories)

	ordered_names = list(model.sanitized_species_names().keys())

	if avg and num_trajectories > 1:
		trajectories = results[0][:,0], np.average(np.dstack(results)[:,1:], axis = 2), ordered_names
	elif not avg or num_trajectories == 1:
		trajectories = results[0][:,0], results[0][:,1:], ordered_names

	#Gillespy2 (since version 1.3) changes the order of the species in the results array to be sorted by alphabet
	#we need to get the indices of rna, prot and the index range for all bound species
	#list(model.sanitized_species_names().keys()).index('rna')

	prot_ind = ordered_names.index('prot')
	rna_ind = ordered_names.index('rna')
	bound_low_ind = ordered_names.index('_' + model.species_names[1])
	bound_high_ind = ordered_names.index('_' + model.species_names[-1])


	#print Kd value based on concentrations at the end of the simulation
	print('Total Kd based on the result from the simulation: ', ((np.mean(pop_to_conc(trajectories[1][-20:-1,prot_ind], volume)) * np.mean(pop_to_conc(trajectories[1][-20:-1,rna_ind], volume))) / (np.mean(pop_to_conc(np.sum(trajectories[1][-20:-1,bound_low_ind:bound_high_ind+1], axis=1), volume)))))


	if plot:
		plot_trajectories(*trajectories)
	
	return trajectories



def plot_trajectories(time, species, names):
	"""
	Creates overview plots from the simulated data.
	INPUT
		trajectories - np.array - results from the simulation
		species - np.array - results from simulation
		names - list, str - species names
	"""

	for i in range(species.shape[1]):
		if names[i] != 'prot':
			plt.plot(time, species[:,i], label=names[i].replace('_',''))


	plt.xlabel('Time')
	plt.ylabel('No. of Species')
	plt.legend()
	plt.tight_layout()
	plt.show()



def conc_to_pop(conc_value, volume):
	"""
	Converts concentrations in Molar into absolute population values
	"""
	return int(conc_value * constants.N_A * volume)



def pop_to_conc(pop_value, volume):
	"""
	Converts absolute population values to concentration in Molar
	"""
	return pop_value/(constants.N_A * volume)



if __name__ == '__main__':
	init_run_print_model('../examples/imp3_1234.csv', num_trajectories = 10)

	#params = get_model_parameters('../examples/imp3_1234.csv')
	#volume = params[5]
	#model = nxn(*params)
	#analytical_kd_result = model.analytical_kd()
	#print('Total Kd based on the result from analytical calculations: ', analytical_kd_result)


