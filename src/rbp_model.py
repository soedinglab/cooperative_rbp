"""
Gillespie model for cooperative binding of a Protein to ssRNA. n binding sites on protein and RNA
"""

# Imports
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
import sys

import gillespy


class nxn(gillespy.Model):
	"""
	Class for defining model parameters, species and reactions. Inherits from the model class of the gillespy library.
	"""

	def __init__(self, (n, prot0, rna0, on, off, volume, lp, d, L)):
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


		self.timespan(np.linspace(0,500,200))


	def gauss_chain(self, L_arg, d_arg):
		#return ((1 + 4*(self.lp/L_arg) + (20/3) * (self.lp/L_arg)**2)/((1-(d_arg/L_arg)**2)**(9/2))*(3/(4*constants.pi *(self.lp/L_arg))**(3/2)) * np.exp(-(3*(d_arg/L_arg)**2)/(4*(self.lp/L_arg)*(1-(d_arg/L_arg)**2))))
		#return ((3/(4*constants.pi*(self.lp/L_arg)))**(3/2)*np.exp(-(3*(d_arg/L_arg)**2)/(4*(self.lp/L_arg)))) #radial distribution function (Becker, Rosa, Everaers, (2010))
		return ((3/(4*constants.pi*self.lp*L_arg))**(3/2)*np.exp(-(3*d_arg**2)/(4*self.lp*L_arg))*self.my_volume*10**(-3)) #radial distribution function of a worm like chain


	def create_reactions(self):
		"""
		Permute all possible reactions and determine the type of reaction, return a list of reactions as gillespy reaction-objects
		"""

		reactions = []

		#loop all species
		for s in self.species_names: 
			
			#find empty binding sites per species
			pos_reac = []
			for ind, elem in enumerate(s):
				if elem == '0':
					pos_reac.append(ind)

			#loop all binding sites per species and determine the possible reaction plus reverse reaction
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
						self.add_parameter(gillespy.Parameter(
							name="chain_distr"+ str(len(reactions)+1),
							expression=self.gauss_chain(L_tot, d_tot)))

						#forward reaction, propensity function (example for 01 -> 11): [01]*[01]*on_stoch_1*chain_distr
						reactions.append(gillespy.Reaction(
							name='r' + str(len(reactions)+1),
							reactants={reactant_object:1},
							products={product_object:1},
							propensity_function= '_' + str(s) + '*' + 'chain_distr' + str(len(reactions)+1) + '*' + 'on_stoch_bi' + str(reac_id) + '*' + '_' + str(s)))

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
						self.add_parameter(gillespy.Parameter(
							name="l_chain_distr"+ str(len(reactions)+1),
							expression=self.gauss_chain(l_L_tot, l_d_tot)))

						self.add_parameter(gillespy.Parameter(
							name="r_chain_distr"+ str(len(reactions)+1),
							expression=self.gauss_chain(r_L_tot, r_d_tot)))

						#forward reaction, propensity function (example for 101 -> 111): [01]*[01]*on_stoch_1*left_chain_distr*right_chain_distr, !!still needs the rates and normalizing!!
						reactions.append(gillespy.Reaction(
							name='r' + str(len(reactions)+1),
							reactants={reactant_object:1},
							products={product_object:1},
							propensity_function= '_' + s + '*' + 'l_chain_distr' + str(len(reactions)+1) + '*' + 'r_chain_distr' + str(len(reactions)+1)))

						#reverse reaction
						reactions.append(gillespy.Reaction(
							name='r' + str(len(reactions)+1),
							reactants={product_object:1},
							products={reactant_object:1},
							rate=self.off_stoch[reac_id]))


		return reactions




def init_run_model(labels=False, num_trajectories=1):
	"""
	Creates an instance of the model, runs the simulation and returns the results.
	INPUT
		labels - bool - turn labels in trajectories on/off
		num_trajectories - int - number of simulations to be returned
	OUTPUT
		list - (trajectories with time and species counts, species names)
	"""
	# Define initial values: 
	# ---------------------
	# n: number of bindings sites on protein and RNA
	# Concentrations [mol L^-1] of unbound Protein and RNA
	# On/Off rate constants for uncooperative binding events (Molar)
	# Volume
	# lp: Persistence length RNA
	# d: Euklidian distance between Protein-Domains (list)
	# L: Distance along RNA chain between the binding sites (list)

	n = 2
	prot0 = 0.2e-6
	rna0 = 0.4e-9
	#off = [1,1]
	on = [3e4, 1.4e5]
	#on = [1.5e6, 9.3e7]
	off = [0.046, 0.13]
	volume = 5.65e-11
	lp = 0.9e-9
	d = np.array([[0, 24.068e-10], [24.068e-10, 0]])
	L = [23*6e-10]#, 23*6.76e-10, 23*6.76e-10]
	model = nxn((n, prot0, rna0, on, off, volume, lp, d, L))
	return (model.run(show_labels=labels, number_of_trajectories=num_trajectories), model.species_names)



def plot_data(time, species, names):
	"""
	Creates plots from the simulated data.
	INPUT
		trajectories - np.array - results from the simulation
	"""
	names.insert(0, 'rna')
	for i in range(species.shape[1]-1):
		plt.plot(time, species[:,i+1], label=names[i])


	plt.xlabel('Time [s]')
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
	results = init_run_model(False, 25)
	time = results[0][0][:,0]
	species = np.dstack(results[0])[:,1:]
	species_avg = np.average(species, axis = 2)

	names = results[1][1:]

	#print((pop_to_conc(species_avg[-1,0]) * pop_to_conc(species_avg[-1,1])) / pop_to_conc(species_avg[-1,-1]))

	plot_data(time, species_avg, names)

	#model = nxn(n, prot0, rna0, on, off, volume, lp, d, L)
	#print(model.serialize())





