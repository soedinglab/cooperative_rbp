"""
Gillespie model for cooperative binding of a Protein to ssRNA. N binding sites on protein and RNA
Includes a class for Gillespie simulations and a method for calculating the Kd based on the analytical solution.
"""

# Imports
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
import sys

import gillespy

#persistence length
#Rechendorff et. al., 2009
lp = 4.6e-9

#length per Base
#Rechendorff et. al., 2009
lpb = 4.3e-10

class nxn(gillespy.Model):
	"""
	Class for defining model parameters, species and reactions. Inherits from the model class of the gillespy library.
	Methods (for external use):
		analytical_kd - calculate the Kd based on the analytical solution
	"""

	def __init__(self, (n, prot0, rna0, on, off, volume, lp, d, L, time)):
		"""
		Initialize the model. Initial values are passed when creating a class instance:
			n: number of bindings sites on protein and RNA
			prot0, rna0: Concentrations [mol L^-1] of unbound Protein and RNA
			On/Off rate constants for uncooperative binding events (Molar) (list)
			Volume
			lp: Persistence length RNA
			d: Euklidian distance between Protein-Domains (np.array, dimension: n,n)
			L: Distance along RNA chain between the binding sites (list)
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


	def gauss_chain(self, x, mu, sig_sq):
		return ((1/(2*constants.pi*sig_sq))**(3/2)*np.exp(-((x-mu)**2)/(2*sig_sq))) #radial distribution function of a worm like chain


	def create_reactions(self):
		"""
		Permute all possible reactions and determine the type of reaction, return a list of reactions as gillespy reaction-objects
		"""

		reactions = []

		#loop all species
		for s in self.species_names: 
			
			#find empty binding sites per species
			pos_reac = [] #list of all empty binding sites in current s
			for ind, elem in enumerate(s):
				if elem == '0':
					pos_reac.append(ind)

			#loop all empty binding sites per species and determine the possible reaction plus reverse reaction
			for pos_reac_ind, reac_id in enumerate(pos_reac):
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

						#total distance between binding sites along RNA chain (L_tot) and between binding domains of the protein (d_tot)
						d_tot = np.nan
						L_tot = np.nan
						if bound_neighbour < reac_id:
							d_tot = self.d[bound_neighbour, reac_id]
							L_tot = sum(self.L[bound_neighbour:(reac_id)])
						elif bound_neighbour > reac_id:
							d_tot = self.d[bound_neighbour, reac_id]
							L_tot = sum(self.L[reac_id:(bound_neighbour)])


						#add gaussian chain distr as parameter
						sig_sq = (2/3) * self.lp * L_tot
						self.add_parameter(gillespy.Parameter(
							name="chain_distr"+ str(len(reactions)+1),
							expression=((self.gauss_chain(d_tot, 0, sig_sq)*10**(-3))/constants.N_A)))

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

						#total distance between binding sites along RNA chain (L_tot) and between binding domains of the protein (d_tot)
						l_d_tot = np.nan
						r_d_tot = np.nan
						l_L_tot = np.nan
						r_L_tot = np.nan

						l_d_tot = self.d[l_bound_neighbour, reac_id]
						l_L_tot = sum(self.L[l_bound_neighbour:reac_id])

						r_d_tot = self.d[r_bound_neighbour, reac_id]
						r_L_tot = sum(self.L[reac_id:r_bound_neighbour])

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
		# loop all species and determine kd for theoretical one step reaction from the unbound state to that species, total kd is the sum of all these indivudal kds
		for s in self.species_names:
			bound_sites = []
			L_tot = []
			d = []
			c_eff = []
			prev_ind = np.nan
			for ind, elem in enumerate(s):
				if elem == '1':
					bound_sites.append(ind)
					if not np.isnan(prev_ind):
						L_tot = sum(self.L[prev_ind:ind])
						d = self.d[prev_ind, ind]
						sig_sq = (2/3) * self.lp * L_tot
						c_eff.append(self.gauss_chain(d, 0, sig_sq) * 10**(-3) / constants.N_A)
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
		params.append(lp) #lp
		params.append(np.fromstring(params_dict['d'], dtype = float, sep = ',').reshape(params[0],params[0])) #d
		params.append([(int(i)*lpb) for i in params_dict['L'].split(',')]) #L
		params.append([int(i) for i in params_dict['time'].split(',')]) #end time, timesteps

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
	
	model = nxn(params)
	print(model.analytical_kd())
	results = model.run(show_labels=labels, number_of_trajectories=num_trajectories)
	if avg and num_trajectories > 1:
		return (results[0][:,0], np.average(np.dstack(results)[:,1:], axis = 2), model.species_names)
	elif not avg or num_trajectories == 1:
		return (results[0][:,0], results[0][:,1:], model.species_names)



def plot_trajectories((time, species, names)):
	"""
	Creates plots from the simulated data.
	INPUT
		trajectories - np.array - results from the simulation
		species - np.array - results from simulation
		names - list, str - species names
	"""
	names.insert(0, 'rna')
	for i in range(species.shape[1]-1):
		plt.plot(time, species[:,i+1], label=names[i+1])


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
	#params = get_model_parameters('../examples/hnrnp_a1.csv')
	#params = get_model_parameters('../examples/zbp1.csv')
	#params = get_model_parameters('../examples/tdp-43.csv')
	params = get_model_parameters('../examples/N_3.csv')
	volume = params[5]


	trajectories = init_run_model(params, num_trajectories = 25)

	#print Kd value based on concentrations at the end of the simulation
	print((pop_to_conc(trajectories[1][-1,0], volume) * pop_to_conc(trajectories[1][-1,1], volume)) / (pop_to_conc(sum(trajectories[1][-1,2:]), volume)))


	plot_trajectories(trajectories)


	#model = nxn(params)
	#print(model.serialize())

