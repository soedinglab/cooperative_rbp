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

import matplotlib.colors as colors

import numpy as np

import rbp_model


def compare_N():
	"""Creates a plot to compare the total Kd for different N, using arbitrary distances and binding constants."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	N = 5
	parameters = 3
	kd = []
	kd_n1 = [1e-5, 1e-5, 1e-5]
	for j in range(1, parameters + 1):
		kd.append(kd_n1[j-1])
		for i in range(2, N+1):
			params = rbp_model.get_model_parameters('../examples/compare_N/' + str(j)  + '/N_' + str(i) + '.csv')
			model = rbp_model.nxn(*params)
			kd.append(model.analytical_kd())

	plot_colors = ['C1', 'b', 'r']
	plot_markers = ['^', 'o', 's']
	plot_labels = ['10nt', '20nt', '40nt']

	for j in range(1, parameters + 1):
		ax.plot(range(1, N+1), kd[((j-1)*N):(j*N)], linestyle='', marker=plot_markers[j-1], color=plot_colors[j-1], label=plot_labels[j-1])

	ax.set_yscale('log')
	ax.set_xticks(range(1,6))
	ax.set_ylabel(r'$K_\text{d}$ [\si{M}]')
	ax.set_xlabel(r'$N$')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(title='RNA linker distance')
	fig.tight_layout()
	fig.savefig('../fig/compare_N.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


def k3_tot(k1, k2, k3):
	params = rbp_model.get_model_parameters('../examples/compare_N//2/N_3.csv')
	model = rbp_model.nxn(*params)
	params[3][0] = 1/k1
	params[3][1] = 1/k2
	params[3][2] = 1/k3
	return model.analytical_kd()
	#call rbp_model.nxn.analytical but not with a parameter file

def compare_kd():
	"""Creates a plot to compare the total Kd for different kd of individual sites, using arbitrary distances and binding constants."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	N = 6
	parameters = 3
	kd = []
	for j in range(1, parameters + 1):
		for i in range(1, N+1):
			params = rbp_model.get_model_parameters('../examples/compare_kd/' + str(j)  + '/N_3_' + str(i) + '.csv')
			model = rbp_model.nxn(*params)
			kd.append(model.analytical_kd())

	plot_colors = ['C1', 'b', 'r']
	plot_markers = ['^', 'o', 's']
	plot_labels = ['\SI{1e-5}{M}', '\SI{1e-4}{M}', '\SI{1e-3}{M}']

	k2 = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 5e-1]

	for j in range(1, parameters + 1):
		ax.plot(k2, kd[((j-1)*N):(j*N)], linestyle='', marker=plot_markers[j-1], color=plot_colors[j-1], label=plot_labels[j-1])

	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_ylabel(r'$K_\text{d}$ [\si{M}]')
	ax.set_xlabel(r'$K_2$ [\si{M}]')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(title='$K_3$:')
	fig.tight_layout()
	fig.savefig('../fig/compare_kd.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


def compare_kd_heat():
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	k1 = 1e-5
	x = np.logspace(-5, -1, 100)
	y = np.logspace(-5, -1, 100)
	z = np.zeros((x.shape[0], y.shape[0]))
	for i, elem_i in enumerate(x):
		for j, elem_j in enumerate(y):
			z[i][j] = k3_tot(k1, elem_j, elem_i)

	xg, yg = np.meshgrid(x,y)
	heat_plot = ax.pcolormesh(x,y,z, norm = colors.LogNorm(vmin=z.min(), vmax=z.max()), cmap = 'jet_r', shading = 'gouraud')

	contour_lines = ax.contour(x,y,z, locator=mpl.ticker.LogLocator(subs=(1,)), linestyles='dashed', linewidths=0.5, colors='black')
	fmt = mpl.ticker.LogFormatterMathtext()
	fmt.create_dummy_axis()
	ax.clabel(contour_lines, inline=False, inline_spacing=0, fmt=fmt, fontsize='smaller')

	clb = fig.colorbar(heat_plot)
	clb.set_label(r'$K_\text{d, tot}(3)$ [\si{M}]')
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_ylabel(r'$K_{\mathrm{d}2}$ [\si{M}]')
	ax.set_xlabel(r'$K_{\mathrm{d}3}$ [\si{M}]')
	#ax.spines['top'].set_visible(False)
	#ax.spines['right'].set_visible(False)
	ax.set_aspect('equal')
	fig.tight_layout()
	fig.savefig('../fig/compare_kd_heat.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


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
	#fig.savefig('../fig/example_overview_no_error.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()



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

	species_names = np.array(['Homo sapiens', 'Mus musculus', 'Drosophila melanogaster', 'Caenorhabditis elegans'])
	no_bins = 7
	#new array size max_value
	domain_count = np.zeros((species_names.shape[0], no_bins+1)) #+1 because of bin with 0 domains

	counter = 0
	# count number of domains, only if species == 'Homo sapiens'
	for ind_i, elem_i in enumerate(domains):
		#if species[ind] == 'Homo sapiens':
		for ind_j, elem_j in enumerate(species_names):
			if species[ind_i] == elem_j:
				counter += 1
				if elem_i < no_bins:
					domain_count[ind_j, elem_i] += 1
				else:
					domain_count[ind_j, -1] += 1


	#plot data
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	for ind, elem in enumerate(domain_count):
		ax.bar(range(1, len(elem)), elem[1:], bottom = np.sum(domain_count[:ind, 1:], axis=0), width = 0.4, label=species_names[ind])

	ax.set_xticks(range(1, domain_count.shape[1]))
	ax.set_xticklabels([1, 2, 3, 4, 5, 6, '7+'])
	ax.set_ylabel(r'RBPs')
	ax.set_xlabel(r'Domains')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend()
	fig.tight_layout()
	fig.savefig('../fig/rbp_distribution.pdf', bbox_inches = 'tight', dpi = 600)
	#plt.show()


if __name__ == '__main__':
	#compare_N()
	#compare_kd()
	compare_kd_heat()
	#N_4_trajectory()
	#N_4_trajectory_detail()
	#example_overview()
	#rbd_distribution()
	pass
