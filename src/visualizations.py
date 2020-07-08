"""
Script to visualize the results from the simulations and anytical calculations.
"""

#Imports

import matplotlib as mpl
#mpl.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern} \usepackage{siunitx}']
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Latin Modern Roman']
#plt.rcParams['pgf.rcfonts'] = False
#plt.rcParams['pgf.preamble'] = ['\\usepackage{lmodern} \\usepackage{siunitx}']

import numpy as np

import rbp_model


def compare_N():
	"""Creates a barplot to compare the total Kd for different N, using arbitrary distances and binding constants."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	kd = []
	for i in range(2, 6):
		params = rbp_model.get_model_parameters('../examples/N_' + str(i) + '.csv')
		model = rbp_model.nxn(*params)
		kd.append(model.analytical_kd())

	kd.insert(0, 1e-5)

	ax.plot(range(1, 6), kd, linestyle='', marker='o', color='black')

	ax.set_yscale('log')
	ax.set_xticks(range(1,6))
	ax.set_ylabel(r'$K_\text{d}$ [\si{M}]')
	ax.set_xlabel(r'$N$')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	fig.tight_layout()
	fig.savefig('../fig/compare_N.pdf', bbox_inches = 'tight', dpi = 600)
	#plt.show()



def N_4_trajectory():
	"""Creates a plot of the trajectory after the simulation with 4 binding sites."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	params = rbp_model.get_model_parameters('../examples/N_3.csv')
	trajectories = rbp_model.init_run_print_model('../examples/N_3.csv', num_trajectories = 5)
	trajectories[2].insert(0, 'rna')
	for i in range(trajectories[1].shape[1]-1):
		ax.plot(trajectories[0], trajectories[1][:,i+1], label = trajectories[2][i+1], linewidth = 1)

	ax.set_ylabel(r'No. of molecules')
	ax.set_xlabel(r'Time')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(ncol = 2)
	fig.tight_layout()
	fig.savefig('../fig/N_3_trajectory.pdf', bbox_inches = 'tight', dpi = 600)
	#plt.show()



def N_4_trajectory_detail():
	"""Creates a zoomed in plot of the trajectory after the simulation with 4 binding sites."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	params = rbp_model.get_model_parameters('../examples/N_4.csv')
	trajectories = rbp_model.init_run_model(params, num_trajectories = 5)
	for i in range(trajectories[1].shape[1]-3):
		ax.plot(trajectories[0], trajectories[1][:,i+2], label = trajectories[2][i+1], linewidth = 1)

	ax.set_ylabel(r'No. of molecules')
	ax.set_xlabel(r'Time')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(ncol = 2)
	fig.tight_layout()
	fig.savefig('../fig/N_4_trajectory_detail.pdf', bbox_inches = 'tight', dpi = 600)
	#plt.show()



def example_overview():
	"""Comparing the calculated/simulated results of the examples to the experimental data."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	#Data
	individual_kd_1 = [2.0e-6, 20.4e-6, 2.1e-6, 9e-6]
	individual_kd_1_error = [0.4e-6, 1.06e-6, 1.3e-6, 2e-6]

	individual_kd_2 = [1.1e-6, 6.8e-6, 2e-6, 4e-6]
	individual_kd_2_error = [0.13e-6, 0.14e-6, 1.3e-6, 2e-6]

	exp_total_kd = [13e-9, 15.5e-9, 10e-9, 33.4e-9]
	exp_total_kd_error = [1e-9, 0.0034e-6, 3e-9, 0.9e-9]
	theoretical_total_kd = [2.4e-9, 17e-9, 7.1e-9, 1.3e-8]
	theoretical_total_kd_error = [1e-10, 6e-10, 9e-10, 0]

	example_count = len(theoretical_total_kd)

	default_style = {"markersize":5, "linestyle":'', "barsabove":True, "ecolor":'black', "capsize":2, "elinewidth":1.5}

	#ax.errorbar(range(1, example_count+1), individual_kd_1, yerr=individual_kd_1_error, marker='o', color='C1', label = 'Individual domains (experimental)', **default_style)
	#ax.errorbar(range(1, example_count+1), individual_kd_2, individual_kd_2_error, marker='o', color='C1', **default_style)
	#ax.errorbar(range(1, example_count+1), exp_total_kd, exp_total_kd_error, marker='^', color='b', label = r'Total $K_\text{d}$ (experimental)', **default_style)
	#ax.errorbar(range(1, example_count+1), theoretical_total_kd, theoretical_total_kd_error, marker='s', color='r', label = r'Total $K_\text{d}$ (calculated)', **default_style)

	ax.plot(range(1, example_count+1), individual_kd_1, marker='o', color='C1', label = 'Individual domains (experimental)', linestyle='')
	ax.plot(range(1, example_count+1), individual_kd_2, marker='o', color='C1', linestyle='')
	ax.plot(range(1, example_count+1), exp_total_kd, marker='^', color='b', label = r'Total $K_\text{d}$ (experimental)', linestyle='')
	ax.plot(range(1, example_count+1), theoretical_total_kd, marker='s', color='r', label = r'Total $K_\text{d}$ (calculated)', linestyle='')
	





	ax.set_ylim(1e-10, 1e-3)
	ax.set_yscale('log')
	ax.set_ylabel(r'$K_\text{d}$ [\si{M}]')
	#ax.set_xlabel()
	ax.set_xticks(range(1,example_count+1))
	ax.set_xticklabels(['ZBP1', 'hnRNP A1', 'PTB34', 'IMP3 - RRM12, KH12'])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(fontsize = 'x-small', loc = 'best')
	fig.tight_layout()
	fig.savefig('../fig/example_overview_no_error.pdf', bbox_inches = 'tight', dpi = 600)



def rbd_distribution():
	"""Showing the distribution of RNA binding domains among RNA binding proteins. Data from RBPDB (http:// rbpdb.ccbr.utoronto.ca/)"""
	def count_domains(domain_string):
		domain_string = str(domain_string)
		if domain_string == 'b\'\\\\N\'' or domain_string == 'b\'\'':
			return 0
		counting = False
		counter = 0
		counter_string = ''
		for elem in domain_string:
			if elem == 'x':
				counting = True
			elif counting == True and elem != ';' and elem != '\'':
				counter_string += elem
			elif (elem == ';' or elem == '\'') and counting == True:
				counter += int(counter_string)
				counter_string = ''
				counting = False

		return counter


	#read data from file
	domains = np.loadtxt('../examples/RBPDB_v1.3.1_proteins.tdt', dtype = int, delimiter = '\t', converters = {8: count_domains}, usecols = (8))
	species = np.loadtxt('../examples/RBPDB_v1.3.1_proteins.tdt', dtype = str, delimiter = '\t',  usecols = (6))

	#find max value
	max_value = 0
	for i in domains:
		if i > max_value:
			max_value = i

	no_bins = 7
	#new array size max_value
	domain_count = np.zeros(no_bins+1) #+1 because of bin with 0 domains

	counter = 0
	# count number of domains, only if species == 'Homo sapiens'
	for ind, elem in enumerate(domains):
		#if species[ind] == 'Homo sapiens':
		counter += 1
		if elem < 7:
			domain_count[elem] += 1
		else:
			domain_count[-1] += 1


	#plot data
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	ax.bar(range(1, len(domain_count)), domain_count[1:], width = 0.4, color='black')

	ax.set_xticks(range(1, len(domain_count)))
	ax.set_xticklabels([1, 2, 3, 4, 5, 6, '7+'])
	ax.set_ylabel(r'RBPs')
	ax.set_xlabel(r'Domains')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	fig.tight_layout()
	fig.savefig('../fig/rbp_distribution.pdf', bbox_inches = 'tight', dpi = 600)
	#plt.show()


if __name__ == '__main__':
	#compare_N()
	#N_4_trajectory()
	#N_4_trajectory_detail()
	example_overview()
	#rbd_distribution()
	pass
