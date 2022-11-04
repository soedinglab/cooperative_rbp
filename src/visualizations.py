"""
Script to visualize the results from the simulations and anytical calculations.
"""

#Imports

import matplotlib as mpl
#mpl.use("pgf")
import matplotlib.pyplot as plt
#mpl.use('AGG')
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{lmodern} \usepackage{siunitx}']
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Latin Modern Roman']
plt.rcParams['mathtext.fontset'] = 'cm'
#plt.rcParams['pgf.rcfonts'] = False
#plt.rcParams['pgf.preamble'] = ['\\usepackage{lmodern} \\usepackage{siunitx}']
plt.rcParams['lines.markersize'] = 4

import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker

import numpy as np
import scipy.optimize
from scipy import constants

from fractions import Fraction

import urllib.parse
import urllib.request

import rbp_model

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

def kd_N():
	"""Plots the total Kd for different N, using arbitrary distances and binding constants."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	N = 5
	parameters = 3
	kd = []
	kd_n1 = [1e-5, 1e-5, 1e-5]
	lengths=[[[10], [10,0,10], [10,0,10,0,10,0,10], [10,0,10,0,10,0,10,0,10]],
		[[20], [20,0,20], [20,0,20,0,20,0,20], [20,0,20,0,20,0,20,0,20]],
		[[40], [40,0,40], [40,0,40,0,40,0,40], [40,0,40,0,40,0,40,0,40]]]
	c_12 = []

	for j in range(1, parameters + 1):
		kd.append(kd_n1[j-1])
		for i in range(2, N+1):
			#params = rbp_model.get_model_parameters('../examples/compare_N/' + str(j)  + '/N_' + str(i) + '.csv')
			#model = rbp_model.nxn(*params)
			kd.append(total_kd(i, L=lengths[j-1][i-2]))
			if i == 2:
				c_12.append(c_eff(i, L=lengths[j-1][i-2]))
			#kd.append(model.analytical_kd())


	plot_colors = plt.cm.Blues(np.linspace(0.5,1,3))
	#plot_colors = ['C1', 'b', 'r']
	plot_markers = ['^', 'o', 's']
	plot_labels = [r'10 nt ($c_{12} = \SI{' + '{:.2e}'.format(c_12[0]) + '}{M}$)', r'20 nt ($c_{12} = \SI{' + '{:.2e}'.format(c_12[1]) + '}{M}$)', r'40 nt ($c_{12} = \SI{' + '{:.2e}'.format(c_12[2]) + '}{M}$)']

	for j in range(1, parameters + 1):
		ax.plot(range(1, N+1), kd[((j-1)*N):(j*N)], linestyle='-', marker=plot_markers[j-1], color=plot_colors[j-1], label=plot_labels[j-1])

	ax.set_yscale('log')
	ax.set_xticks(range(1,6))
	ax.set_ylabel(r'${K_\mathrm{av}(n)}^{-1}$ [M]')
	ax.set_xlabel(r'$n$')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(title='RNA linker distance')
	fig.tight_layout()
	fig.savefig('../fig/kd_N.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


def kd_tot_k3():
	"""Plots the total Kd for different kd of individual sites, using arbitrary distances and binding constants."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	N = 6
	parameters = 3
	k1=1e-5
	k3 = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 5e-1]
	k2 = [1e-5, 1e-4, 1e-3]
	kd = []
	for j in range(1, parameters + 1):
		for i in range(1, N+1):
			pass
			#params = rbp_model.get_model_parameters('../examples/compare_kd/' + str(j)  + '/N_3_' + str(i) + '.csv')
			#model = rbp_model.nxn(*params)
			#kd.append(model.analytical_kd())


	for i in k2:
		for j in k3:
			kd.append(total_kd(3, kd=np.array([k1,i,j])))

	plot_colors = plt.cm.Greens(np.linspace(0.5,1,3))
	#plot_colors = ['C1', 'b', 'r']
	plot_markers = ['^', 'o', 's']
	plot_labels = ['$10^{-5}$ M', '$10^{-4}$ M', '$10^{-3}$ M']


	for j in range(1, parameters + 1):
		ax.plot(k3, kd[((j-1)*N):(j*N)], linestyle='-', marker=plot_markers[j-1], color=plot_colors[j-1], label=plot_labels[j-1])


	params = rbp_model.get_model_parameters('../examples/N_3.csv')
	model = rbp_model.nxn(*params)
	c_eff = model.get_concentration(1,2)

	ax.vlines(c_eff, 1e-9, 1e-5, ls ='dashed', linewidth=1)
	ax.text(c_eff- 2e-4, 1.2e-5, r'$c_{23}$')

	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_ylabel(r'${K_\mathrm{av}}^{-1}$ [M]')
	ax.set_xlabel(r'$K_{\mathrm{d,}3}$ [M]')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(title='$K_{\mathrm{d,}2}$:')
	fig.tight_layout()
	fig.savefig('../fig/kd_tot_k3.pdf', bbox_inches = 'tight', dpi = 600)
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
	fig, ax = plt.subplots(1,1, figsize=(4.6,3))

	#Data
	individual_kd_1 = [2.0e-6, 20.4e-6, 2.1e-6, 9e-6, 1.75e-6, 5e-3, 5e-3, 5e-3, 390e-6]

	individual_kd_2 = [1.1e-6, 6.8e-6, 2e-6, 4e-6, 1.4e-4, 80e-6, 80e-6, 80e-6, 140e-6]
	individual_kd_3 = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 350e-6]

	exp_total_kd = np.array([13e-9, 15.5e-9, 10e-9, 33.4e-9, 4.95e-8, 16.4e-6, 7.1e-6, 1.32e-6, 7.3e-9])
	theoretical_total_kd = [2.4e-9, 17e-9, 7.1e-9, 1.3e-8, 1.1e-7, 4.1e-5, 2.2e-5, 3.6e-6, 27e-9]

	indep_binding_tot_kd = (np.array(individual_kd_1)**(-1) + np.array(individual_kd_2)**(-1))**(-1)

	example_count = len(theoretical_total_kd)

	default_style = {"markersize":5, "linestyle":'', "barsabove":True, "ecolor":'black', "capsize":2, "elinewidth":1.5}

	x_values = np.array([1, 2, 3, 4, 5, 6, 6.5, 7, 8])
	x_offset = 0.05
	ax.plot(x_values+[0,0,0,0,x_offset,x_offset,x_offset,x_offset,0], indep_binding_tot_kd, marker='^', color=plt.cm.tab20([2])[0], label = r'$K_\mathrm{d}$ indep. binding (calculated)', linestyle='')
	ax.plot(x_values+[0,0,-x_offset,0,-x_offset,0,0,0,-x_offset], individual_kd_1, marker='o', color=plt.cm.tab20([18])[0], label = 'Individual domains (experimental)', linestyle='')
	ax.plot(x_values+[0,0,x_offset,0,0,-x_offset,-x_offset,-x_offset,0], individual_kd_2, marker='o', color=plt.cm.tab20([18])[0], linestyle='')
	ax.plot(x_values+[0,0,0,0,0,0,0,0,x_offset], individual_kd_3, marker='o', color=plt.cm.tab20([18])[0], linestyle='')
	ax.plot(x_values+[0,-x_offset,0,0,0,0,0,0,0], exp_total_kd, marker='D', color=plt.cm.tab20([0])[0], label = r'Total $K_\mathrm{d}$ (experimental)', linestyle='')
	ax.plot(x_values+[0,x_offset,0,0,0,0,0,0,0], theoretical_total_kd, marker='s', color=plt.cm.tab20([6])[0], label = r'Total $K_\mathrm{d}$ (calculated)', linestyle='')
	
	annotation_text = ['8\,nt', '4\,nt', '1\,nt']
	for i, j in enumerate([5, 6, 7]):
		plt.annotate(annotation_text[i], (x_values[j], exp_total_kd[j] / 10), ha='center', va='center')

	ax.set_ylim(1e-10, 1e-2)
	ax.set_yscale('log')
	ax.set_ylabel(r'$K_\mathrm{d}$ [M]')
	#ax.set_xlabel()
	ax.set_xticks([1, 2, 3, 4, 5, 6.5, 8])
	ax.set_xticklabels(['ZBP1', 'hnRNP A1', 'PTB34', r'IMP3\\RRM12\\KH12', r'IMP1\\KH1,KH2', 'U2AF', r'KSRP\\KH1 Mutant'], rotation=45)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(fontsize = 'small', loc = (0.01,0.8))
	fig.tight_layout()
	fig.savefig('../fig/example_overview_revision.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


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


def rbd_distribution():
	"""Showing the distribution of RNA binding domains among RNA binding proteins. Data from RBPDB (http:// rbpdb.ccbr.utoronto.ca/)"""


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
	#new array size no_bins or max_value
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
	plot_colors = plt.cm.Set3([0,5,2,3])

	fig, ax = plt.subplots(1,1, figsize=(4.7,2.82))

	for ind, elem in enumerate(domain_count):
		ax.bar(range(1, len(elem)), elem[1:], bottom = np.sum(domain_count[:ind, 1:], axis=0), width = 0.4, label=species_names[ind], color=plot_colors[ind])

	ax.set_xticks(range(1, domain_count.shape[1]))
	ax.set_xticklabels([1, 2, 3, 4, 5, 6, '7+'])
	ax.set_ylabel(r'RBPs')
	ax.set_xlabel(r'Domains')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend()
	fig.tight_layout()
	fig.savefig('../fig/rbp_distribution.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


def oligomer_distribution():
	"""Showing the distribution of Mono-, Di- and Oligomers among RNA binding proteins. Data from PDBePISA (https://www.ebi.ac.uk/pdbe/pisa/)."""


	domains = np.loadtxt('../examples/RBPDB_v1.3.1_proteins.tdt', dtype = int, delimiter = '\t', converters = {8: count_domains}, usecols = (8))
	uniprot_ids = np.loadtxt('../examples/RBPDB_v1.3.1_proteins.tdt', dtype = str, delimiter='\t', usecols = (14))

	total_proteins = 0
	
	#Map UniProt IDs to PDB Ids
	uniprot_ids_str = ''
	count = 0
	for i, elem in enumerate(domains):
		if elem == 1:
			total_proteins += 1
		if elem == 1 and uniprot_ids[i] != '\\N':
			uniprot_ids_str += ' ' + uniprot_ids[i].replace(';', ' ')
			count += 1

	url = 'https://www.uniprot.org/uploadlists/'

	params = {
	'from': 'ID',
	'to': 'PDB_ID',
	'format': 'tab',
	'query': uniprot_ids_str,
	}

	data = urllib.parse.urlencode(params)
	data = data.encode('utf-8')
	req = urllib.request.Request(url, data)
	with urllib.request.urlopen(req) as f:
		response = f.read()
	#print(response.decode('utf-8').split('\n'))

	id_map = [[],[],[]]
	
	#parse output of mapping of uniprot ids to pdb ids
	for elem in response.decode('utf-8').split('\n')[1:]:
		if elem.split('\t') != [''] and elem.split('\t')[1]:
			id_map[0].append(elem.split('\t')[0])
			id_map[1].append(elem.split('\t')[1])

	print('From ' + str(total_proteins) + ' proteins with one domain, we have ' + str(len(uniprot_ids_str.split())) + ' UniProt Ids, which can be connected to ' + str(len(id_map[0])) + ' PDB Ids.')

	#quataernary structure data
	qs_ids = np.loadtxt('../examples/pdbepisa.dat', dtype =str, usecols=(1), skiprows = 2)
	qs_subunits = np.loadtxt('../examples/pdbepisa.dat', dtype =int, usecols=(2), skiprows = 2)

	#match qs data to the pdbs in our dataset
	for i in range(len(id_map[0])):
		match_ids = np.where(qs_ids==id_map[1][i].lower())
		for j in match_ids[0]:
			id_map[2].append(qs_subunits[j])
		if match_ids[0].size==0:
			id_map[2].append(np.nan)


	#for i in range(0,len(id_map[0])-1):
	#	print(id_map[0][i], id_map[2][i])

	#count number of oligomers, remove redundant PDB Ids
	nans = 0
	non_redundant_pdbs = []
	no_bins = 7
	oligomer_count = np.zeros((no_bins))
	for i, elem in enumerate(id_map[2]):
		#print(id_map[0][i], id_map[1][i], id_map[2][i])
		if id_map[1][i] in non_redundant_pdbs:
			continue
		if elem > (no_bins - 1):
			oligomer_count[-1] += 1
		elif np.isnan(elem):
			nans += 1
			non_redundant_pdbs.append(id_map[1][i])
			continue
		else:
			oligomer_count[elem-1] += 1
		non_redundant_pdbs.append(id_map[1][i])
	
	#count non redundant UniProt Ids in the result structures
	non_redundant_result_uniprots = []
	for elem in non_redundant_pdbs:
		index = [i for i, e in enumerate(id_map[1]) if e == elem]
		for i in index:
			if id_map[0][i] not in non_redundant_result_uniprots:
				non_redundant_result_uniprots.append(id_map[0][i])



	
	print('The UniProt Ids could be connected to ' + str(np.sum(oligomer_count) + nans) + ' unique PDB Ids. Out of these, the PDBePISA database contained data for ' + str(np.sum(oligomer_count)) + ' protein structures, from ' + str(len(non_redundant_result_uniprots)) + ' unique proteins.')

	fig, ax = plt.subplots(1,1, figsize=(4,2.4))
	ax.bar(range(1, no_bins+1), oligomer_count,  width = 0.4, color='black')

	ax.set_xticks(range(1, no_bins+1))
	ax.set_xticklabels([1, 2, 3, 4, 5, 6, '7+'])
	ax.set_ylabel(r'PDB structures')
	ax.set_xlabel(r'Oligomeric state')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	#ax.legend()
	fig.tight_layout()
	fig.savefig('../fig/oligomer_distribution.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


def kd_motif_density():
	"""Plots the dependence of the Kd of the density of binding motifs on the RNA."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	RNA_length = 200
	RNA_conc = 1e-7
	kd_1 = 50e-6

	linker_length = np.zeros(7)
	linker_length[0] = RNA_length
	for i in range(1, linker_length.shape[0]):
		linker_length[i] = int(linker_length[i-1] /2)
	motif_density = 1 / linker_length
	motif_count = np.ceil(motif_density * RNA_length)
	print(linker_length)
	print(motif_count)
	print(motif_density)

	kd_2 = []
	kd_3 = []
	kd_4 = []
	#kd_6 = []
	for i, elem in enumerate(linker_length):
		kd_2.append(total_kd(2, L = np.array([elem]), kd=np.array([kd_1, kd_1])))
		kd_3.append(total_kd(3, L = np.array([elem, 0, elem]), kd=np.array([kd_1, kd_1, kd_1])))
		kd_4.append(total_kd(4, L = np.array([elem, 0, elem, 0, elem]), kd=np.array([kd_1, kd_1, kd_1, kd_1])))
		#kd_6.append(total_kd(6, L = np.array([elem, 0, elem, 0, elem, 0, elem, 0, elem]), kd=np.array([kd_1, kd_1, kd_1, kd_1, kd_1, kd_1])))


	kd_1_density = [kd_1 / i for i in motif_count]
	kd_2_density = [elem / (motif_count[ind] - 1) for ind, elem in enumerate(kd_2[1:], start=1)]
	kd_3_density = [elem / (motif_count[ind] - 2) for ind, elem in enumerate(kd_3[2:], start=2)]
	kd_4_density = [elem / (motif_count[ind] - 3) for ind, elem in enumerate(kd_4[2:], start=2)]
	#kd_6_density = [elem / (motif_count[ind] - 5) for ind, elem in enumerate(kd_6[3:], start=3)]

	# 1 binding site, 2 domains
	kd_2_density.insert(0, (kd_1/2))
	
	# 2 and 1 binding sites, 3 domains
	kd_3_density.insert(0, (kd_2[1]/2))
	kd_3_density.insert(0, (kd_1/3))

	# 2 and 1 binding sites, 4 domains
	kd_4_density.insert(0, (kd_2[1]/3))
	kd_4_density.insert(0, (kd_1/4))

	# 4 and 2 and 1 binding sites, 6 domains 
	#kd_6_density.insert(0, (kd_4[2]/3))
	#kd_6_density.insert(0, (kd_2[1]/5))
	#kd_6_density.insert(0, (kd_1/6))


	colors = plt.cm.Purples(np.linspace(0.5,1,5))
	
	#occupancy_1 = RNA_conc/(1e-5+RNA_conc)
	#ax.hlines(occupancy_1, 0, 100, linestyles = 'dashed', linewidth = 1, label='1')
	ax.plot(motif_density, kd_1_density, linestyle='-', label='   ', color=colors[0], marker='^')
	ax.plot(motif_density, kd_2_density, linestyle='-', label='   ', color=colors[1], marker='o')
	ax.plot(motif_density, kd_3_density, linestyle='-', label='   ', color=colors[2], marker='s')
	ax.plot(motif_density, kd_4_density, linestyle='-', label='   ', color=colors[3], marker='D')
	#ax.plot(motif_density, kd_6_density, linestyle='-', label='6', color=colors[4], marker='h')

	RNA_conc = 0.1e-6
	ax.hlines(RNA_conc, (1/200), (1/3), linestyle = 'dotted', linewidth=1, color='black')
	ax.text(1/200, 8e-8, r'[RNA]$ = \SI{0.1}{\micro M}$', ha='left', va='top')
	#ax.set_yticks([0.1e-6], minor=True)
	#ax.set_yticklabels([r'[RNA]$ = \SI{0.1}{\micro M}$'], minor=True)

	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.xticks(motif_density, [str(Fraction(i).limit_denominator()) for i in motif_density])
	ax.set_xticks([], minor=True)

	ax.set_ylabel(r'${K_\mathrm{av}}^{-1}$ [M]')
	ax.set_xlabel(r'Binding site density [nt$^{-1}$]')
	#ax.set_ylim(-0.1,1.1)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(title='No. of RBDs', loc='best')
	fig.tight_layout()
	fig.savefig('../fig/kd_motif_density.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


def occupancy_motif_density():
	"""Plots the dependence of the RNA occupancy of the density of binding motifs on the RNA."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	RNA_length = 200
	RNA_conc = 1e-7
	kd_1 = 50e-6

	linker_length = np.zeros(7)
	linker_length[0] = RNA_length
	for i in range(1, linker_length.shape[0]):
		linker_length[i] = int(linker_length[i-1] /2)
	motif_density = 1 / linker_length
	motif_count = np.ceil(motif_density * RNA_length)
	print(linker_length)
	print(motif_count)
	print(motif_density)

	kd_2 = []
	kd_3 = []
	kd_4 = []
	#kd_6 = []
	for i, elem in enumerate(linker_length):
		kd_2.append(total_kd(2, L = np.array([elem ]), kd=np.array([kd_1, kd_1])))
		kd_3.append(total_kd(3, L = np.array([elem, 0, elem]), kd=np.array([kd_1, kd_1, kd_1])))
		kd_4.append(total_kd(4, L = np.array([elem, 0, elem, 0, elem]), kd=np.array([kd_1, kd_1, kd_1, kd_1])))
		#kd_6.append(total_kd(6, L = np.array([elem, 0, elem, 0, elem, 0, elem, 0, elem]), kd=np.array([kd_1, kd_1, kd_1, kd_1, kd_1, kd_1])))


	occupancy_1 = [(RNA_conc/(kd_1 / i + RNA_conc)) for i in motif_count]
	occupancy_2 = [(RNA_conc/(elem / (motif_count[ind] - 1) + RNA_conc)) for ind, elem in enumerate(kd_2[1:], start=1)]
	occupancy_3 = [(RNA_conc/(elem / (motif_count[ind] - 2) + RNA_conc)) for ind, elem in enumerate(kd_3[2:], start=2)]
	occupancy_4 = [(RNA_conc/(elem / (motif_count[ind] - 3) + RNA_conc)) for ind, elem in enumerate(kd_4[2:], start=2)]
	#occupancy_6 = [(RNA_conc/(elem / (motif_count[ind] - 5) + RNA_conc)) for ind, elem in enumerate(kd_4[3:], start=3)]

	# 1 binding site, 2 domains
	occupancy_2.insert(0, (RNA_conc /(RNA_conc + kd_1/2)))
	
	# 2 and 1 binding site, 3 domains
	occupancy_3.insert(0, (RNA_conc /(RNA_conc + kd_2[1]/2)))
	occupancy_3.insert(0, (RNA_conc /(RNA_conc + kd_1/3)))

	# 2 and 1 binding sites, 4 domains
	occupancy_4.insert(0, (RNA_conc / (RNA_conc + kd_2[1]/4)))
	occupancy_4.insert(0, (RNA_conc / (RNA_conc + kd_1/4)))

	# 4 and 2 and 1 binding sites, 6 domains
	#occupancy_6.insert(0, (RNA_conc / (RNA_conc + kd_4[2]/3)))
	#occupancy_6.insert(0, (RNA_conc / (RNA_conc + kd_2[1]/5)))
	#occupancy_6.insert(0, (RNA_conc / (RNA_conc + kd_1/6)))

	colors = plt.cm.Purples(np.linspace(0.5,1,5))

	####
	# fit to hill curve
	motif_density_small_steps = np.logspace(np.log2(1/200), np.log2(1/3), num=100, base=2)
	guess = (1,1)
	hill_param_1 = hill_fit(occupancy_1, motif_density, guess)
	print('Parameters for fit to Hill function with 1 domain: ', hill_param_1)
	hill_plot = plt.plot(motif_density_small_steps, hill_func(motif_density_small_steps, *hill_param_1[0]), color=colors[0])

	guess = (0.01,1)
	hill_param_2 = hill_fit(occupancy_2, motif_density, guess)
	print('Parameters for fit to Hill function with 2 domains: ', hill_param_2)
	hill_plot = plt.plot(motif_density_small_steps, hill_func(motif_density_small_steps, *hill_param_2[0]), color=colors[1])

	guess = (0.01,2)
	hill_param_3 = hill_fit(occupancy_3, motif_density, guess)
	print('Parameters for fit to Hill function with 3 domains: ', hill_param_3)
	hill_plot = plt.plot(motif_density_small_steps, hill_func(motif_density_small_steps, *hill_param_3[0]), color=colors[2]) 
	guess = (0.0001,3)
	hill_param_4 = hill_fit(occupancy_4, motif_density, guess)
	print('Parameters for fit to Hill function with 4 domains: ', hill_param_4)
	hill_plot = plt.plot(motif_density_small_steps, hill_func(motif_density_small_steps, *hill_param_4[0]), color=colors[3])

	#guess = (0.0001,3)
	#hill_param_6 = hill_fit(occupancy_6, motif_density, guess)
	#print('Parameters for fit to Hill function with 6 domains: ', hill_param_6)
	#hill_plot = plt.plot(motif_density_small_steps, hill_func(motif_density_small_steps, *hill_param_6[0]), color=colors[4])


	#occupancy_1 = RNA_conc/(1e-5+RNA_conc)
	#ax.hlines(occupancy_1, 0, 100, linestyles = 'dashed', linewidth = 1, label='1')
	ax.plot(motif_density, occupancy_1, linestyle='', label='   ', color=colors[0], marker='^')
	ax.plot(motif_density, occupancy_2, linestyle='', label='   ', color=colors[1], marker='o')
	ax.plot(motif_density, occupancy_3, linestyle='', label='   ', color=colors[2], marker='s')
	ax.plot(motif_density, occupancy_4, linestyle='', label='   ', color=colors[3], marker='D')
	#ax.plot(motif_density, occupancy_6, linestyle='', label='6', color=colors[4], marker='h')

	ax.hlines(0.5, (1/200), (1/3), linestyle='dashdot', linewidth=1)

	ax.set_xscale('log')
	#ax.set_yscale('log')
	plt.xticks(motif_density, [str(Fraction(i).limit_denominator()) for i in motif_density])
	ax.set_xticks([], minor=True)

	ax.set_ylabel(r'Relative occupancy')
	ax.set_xlabel(r'Binding site density [nt$^{-1}$]')
	#ax.set_ylim(-0.1,1.1)
	#ax.set_xlim(1/250, 1/2.5)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.legend(title=r'\begin{center}[RNA]$ = \SI{0.1}{\micro M}$\\No. of RBDs\end{center}', loc=[0.02, 0.52])
	fig.tight_layout()
	fig.savefig('../fig/occupancy_motif_density_fit.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


def kd_lrna_lprot():
	"""Plots the Kav in dependence of RNA and protein linker lengths."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	L = np.arange(1, 20, 1)
	L_p = np.arange(1, 100, 1)
	kd = np.zeros((L.shape[0], L_p.shape[0]))

	for i, elem_i in enumerate(L):
		for j, elem_j in enumerate(L_p):
			kd[i,j] = total_kd(2, L = np.array([elem_i]), L_p = np.array([elem_j]))

	color_plot = plt.pcolormesh(np.append(L_p-0.5, L_p[-1]+0.5), np.append(L-0.5, L[-1]+0.5), kd, norm=colors.LogNorm(vmin=kd.min(), vmax=kd.max()), cmap='magma_r', rasterized=True)
	cbar = plt.colorbar(color_plot, label=r'${K_\mathrm{av}}^{-1}$ [M]')
	
	plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
	plt.gca().set_xlabel(r'$L_\text{protein}$ [aa]')
	plt.gca().set_ylabel(r'$L$ [nt]')
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	fig.tight_layout()
	#fig.savefig('../fig/kd_lrna_lprot.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


def ceff_L():
	"""Plots the effective concentration (particle in a sphere and gaussian/worm-like chain) for different RNA linker lengths between binding motifs."""
	fig, ax = plt.subplots(1,1, figsize=(5,3))

	L = np.arange(1, 35, 1)
	d = np.array([[0, 3e-9], [3e-9, 0]])
	c_gauss = np.zeros((L.shape[0]))
	c_sphere = ((4/3) * constants.pi * L**3)**(-1)
	kd_1 = 1e-5

	for i, elem_i in enumerate(L):
			c_gauss[i] = c_eff(2, L = np.array([elem_i]), d = d, kd=np.array([kd_1,kd_1]))

	
	plt.plot(L, c_gauss*10**(3), label='WLC', color='red', marker='o', linestyle='dashed')
	plt.plot(L[3:], c_sphere[3:]*10**(3), label='Uniform', color='black', marker='s', linestyle='dashed')

	plt.yscale('log')

	#plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))
	plt.gca().set_xlabel(r'$l_{12}$ [nt]')
	plt.gca().set_ylabel(r'$c_{12}$ [mM]')
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.legend()
	fig.tight_layout()
	fig.savefig('../fig/ceff_L.pdf', bbox_inches = 'tight', dpi = 600)
	plt.show()


def total_kd(N, **kwargs):
	"""Return the total Kd for N binding sites with the given parameters."""
	params = rbp_model.get_model_parameters('../examples/N_' + str(N) + '.csv')
	#params: n, prot0, rna0, on (1/kd), off, volume, lp, d, L, L_p, lp_p, time
	#change parameters according to the parameters passed in kwargs
	parameter_map = {'kd': 3, 'd': 7, 'L': 8, 'L_p': 9}
	for key, value in kwargs.items():
		param_index = parameter_map.get(key)
		if key == 'kd':
			params[param_index] = 1/value
		if key == 'd':
			params[param_index] = value
		if key == 'L':
			params[param_index] = [i*lpb for i in value]
			#params[param_index] = value * lpb
		if key == 'L_p':
			params[param_index] = [i*lpaa for i in value]


	model = rbp_model.nxn(*params)
	return model.analytical_kd()

def c_eff(N, **kwargs):
	"""Return the effective local concentration for N binding sites with the given parameters."""
	params = rbp_model.get_model_parameters('../examples/N_' + str(N) + '.csv')
	#change parameters according to the parameters passed in kwargs
	parameter_map = {'kd': 3, 'd': 7, 'L': 8, 'L_p': 9}
	for key, value in kwargs.items():
		param_index = parameter_map.get(key)
		if key == 'kd':
			params[param_index] = 1/value
		if key == 'd':
			params[param_index] = value
		if key == 'L':
			params[param_index] = [i*lpb for i in value]
			#params[param_index] = value * lpb
		if key == 'L_p':
			params[param_index] = [i*lpaa for i in value]


	model = rbp_model.nxn(*params)
	return model.get_concentration(0,1)


def hill_func(conc, ka, n):
	return(conc**n / (ka**n + conc**n))

def hill_fit(occupancy, motif_density, guess):
	return scipy.optimize.curve_fit(hill_func, motif_density, occupancy, p0=(guess[0], guess[1]))

if __name__ == '__main__':
	#kd_N()
	#kd_tot_k3()
	#N_4_trajectory()
	#N_4_trajectory_detail()
	example_overview()
	#rbd_distribution()
	#oligomer_distribution()
	#kd_motif_density()
	#occupancy_motif_density()
	#kd_lrna_lprot()
	#ceff_L()
