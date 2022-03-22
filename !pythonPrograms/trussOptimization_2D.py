#import
import designLibrary_test
import trussAnalysis_2D
import random
import copy
import keyboard
import matplotlib.pyplot as plt
import time
import numpy as np

#functions
def gather_node_info():
	all_node_info = []
	select_node_info = []
	num_elements = len(str_values[0])
	ele_values = str_values[0]
	for ele_props in ele_values:
		n1 = ele_props[0][0]
		n2 = ele_props[0][1]
		if n1 not in all_node_info:
			all_node_info.append(n1)
		if n2 not in all_node_info:
			all_node_info.append(n2)
	#determine which nodes can be altered...
	#any free, unloaded nodes can be altered by GA
	for node_info in all_node_info:
		n1x = node_info[0][0]  ;  n1z = node_info[0][2]
		nsx = node_info[1][0]  ;  nsz = node_info[1][2]
		nlx = node_info[2][0]  ;  nlz = node_info[2][2]
		if nsx == 1 and nsz == 1: #if the node is free,
			if nlx == 0 and nlz == 0: #if the node has no loads,
				select_node_info.append(node_info[0])
	genes = [0] * ( len(select_node_info)*2 + num_elements*2 )
	return select_node_info, genes

	#temporary print to check function. delete when finished.
	for node_info in select_node_info:
		print(node_info)
		
def calculate_fitness(str_values):
	disp = max(str_values[10]) #maximum nodal displacement
	dfit = -(disp/dref)+1
	
	#minimum safety factor
	sf_target = 1.0
	sf_min = min(str_values[14])
	if sf_min < sf_target:
		sfit = (sf_min/sf_target)-.5
	if sf_min >= sf_target:
		sfit = 1-(sf_min - sf_target)/sf_target
		
	# minimize mass
	mass = str_values[13]
	mfit = 1 - (mass / mref)

	# score multipliers
	dfit = dfit * 1
	sfit = sfit * 1
	mfit = mfit * 3.5

	fitness = dfit + sfit + mfit
	return fitness, dfit, sfit, mfit
	
def mutate_genes():	
	changes = []
	for i in range(0, max_mutations):
		mutate = []
		random_gene = random.randint(0, len(genes)-1)
		random_change = random.uniform(-gene_range, gene_range)
		mutate.append(random_gene)
		mutate.append(random_change)
		changes.append(mutate)
		genes[random_gene] += random_change
	return genes, changes

def insert_genes():
	child = copy.deepcopy(str_values_ini)
	ele_values = child[0]
	for ele_props in ele_values:
		#define the coordinates of the two nodes of an element
		node1 = ele_props[0][0][0]
		node2 = ele_props[0][1][0]
		n1x = node1[0]  ;  n1z = node1[2]
		n2x = node2[0]  ;  n2z = node2[2]
		for i in range(0, len(variable_nodes)):
			#for each element, run thru the list of variable nodes
			v_node = variable_nodes[i]
			v_node_x = v_node[0]
			v_node_z = v_node[2]
			gene_x = genes[2*i] # 0, 2, 4, etc.
			gene_z = genes[2*i+1] # 1, 3, 5, etc.
			if v_node == node1:
				#if node1 of the element is variable,
				#mutate with the corresponding genes.
				ele_props[0][0][0][0] += gene_x
				ele_props[0][0][0][2] += gene_z
			if v_node == node2:
				#if node2 of the element is variable,
				#mutate with the corresponding genes.
				ele_props[0][1][0][0] += gene_x
				ele_props[0][1][0][2] += gene_z
	#mutate element xc_area with the mutated genes
	num_genes = len(genes)
	num_elements = len(ele_values)
	num_coord_genes = num_genes - num_elements*2
	for i in range(0, num_elements):
		gene_xcarea = genes[num_coord_genes + i]/1000#reduction factor
		ele_values[i][1] += gene_xcarea
		#prevent very small and negative xc areas
		min_area = 0.3 /10000 #convert sq.cm to sq.m
		if ele_values[i][1] <= min_area:
			ele_values[i][1] = min_area #1 sq.cm
	for i in range(0, num_elements):
		skip_genes = num_coord_genes + num_elements
		gene = genes[skip_genes + i]
		gr = gene_range
		if gene < (-1/3)*gr:
			gene_mat = 'a' #steel
		if (-1/3)*gr <= gene and gene < (1/3)*gr:
			gene_mat = 's' # aluminum
		if (1/3)*gr <= gene:
			gene_mat = 't' # titanium
		ele_values[i][2] = gene_mat
	return child

def undo_mutate_genes():
	for i in range(0, len(changes)):
		mutation = changes[i]
		random_gene = mutation[0]
		random_change = mutation[1]
		genes[random_gene] -= random_change
	return genes
	
def print_results():
	ele_values_p = parent[0]
	print("\nxc area, material code, n1 coords, n2 coords, FOS,fMode")
	for i in range(0, len(ele_values_p)):
		ele_props = ele_values_p[i]
		xc_area = ele_props[1]
		mat_code = ele_props[2]
		fos = parent[14][i] #Factor of Safety
		fm = parent[15][i] #failure mode
		total_mass = parent[13]
		max_disp = max(parent[10]) * 1000 # convert m to mm
		n1_x = ele_props[0][0][0][0]
		n1_z = ele_props[0][0][0][2]
		n2_x = ele_props[0][1][0][0]
		n2_z = ele_props[0][1][0][2]
		print("% .6f" % xc_area, mat_code,"% .3f" % n1_x,
			"% .3f" % n1_z,"% .3f" % n2_x,"% .3f" % n2_z,
			 "% .3f" % fos, fm)
	print('total structure mass: '+str("%.3f" % total_mass)+' kg')
	print('max structural displacement: '+str("%.3f" % max_disp)+' mm')
	
def info_plotter():
	figTitle = '2D Truss Optimization Using a Genetic Algorithm'
	info_fig.suptitle(figTitle, fontsize=18, color='darkblue')
	sub1 = info_fig.add_subplot(2,3,1)
	sub1.plot(x, y_parent, color='blue', label = 'parent fitness')
	sub1.scatter(x, y_child, color='red', s=.5, marker ='.',
	label='child fitness')
	sub1.scatter(x, y_dfit, color='darkcyan', s=.5, marker ='.',
		label='dfit score')
	sub1.scatter(x, y_sfit, color='darkorange', s=.5, marker ='.',
		label='sfit score')
	sub1.scatter(x, y_mfit, color='darkred', s=.5, marker ='.',
		label='mfit score')	
	sub1.set(xlabel = 'iterations', ylabel = 'fitness value')
	sub1.legend(loc='upper center', bbox_to_anchor=(.3,1.3), ncol=2,
		markerscale=20)	
	sub1.grid()
	
	sub2 = info_fig.add_subplot(2,3,2)
	s2color = 'darkgoldenrod'
	sub2.scatter(x, y_mass, color=s2color, s=5, marker ='.')
	sub2.set(xlabel='iterations', ylabel='structure mass, kg')
	sub2.grid()

	sub4 = info_fig.add_subplot(2,3,4)
	s4color = 'green'
	sub4.scatter(x, y_disp, color = s4color, s=5, marker ='.')
	sub4.set(xlabel='iterations', ylabel='max displacement, mm')
	sub4.grid()

	sub5 = info_fig.add_subplot(2,3,5)
	minColor = 'mediumblue'
	sub5.scatter(x, y_fos_min, color=minColor, s=10, marker ='.')
	sub5.set_xlabel('iterations', color = 'black')
	sub5.set_ylabel('MIN Factor of Safety (dot)', color = minColor)
	sub5.tick_params(axis='y', colors = minColor, labelsize='medium')
	sub5.grid()

	sub5a = sub5.twinx()
	maxColor = 'purple'
	sub5a.scatter(x, y_fos_max, color=maxColor, s=10, marker ='x')
	sub5a.set_ylabel('MAX Factor of Safety (x)', color=maxColor)
	sub5a.tick_params(axis = 'y', colors = maxColor)
	
	sub3 = info_fig.add_subplot(2,3,3)
	sub3.axis('off')
	#create a list of z axis spacing values for the text summary
	text_Z_spacing = 0.057 ; z = np.arange(1, -1, -text_Z_spacing)
	sub3.text(0,z[0],'Optimization Summary-',fontweight='bold',
		style='italic',size=14, color='blue')
	sub3.text(0,z[1],'')
	sub3.text(0,z[2],'Original Design-',fontweight='bold')
	sub3.text(0,z[3],'Mass:'+str("% .3f" % orig_mass)+' kg')
	sub3.text(0,z[4],'Max. Disp.:'+str("% .2f" % orig_max_disp)+' mm')
	sub3.text(0,z[5],'Min. FoS:'+str("% .2f" % orig_min_fos))
	sub3.text(0,z[6],'Max. FoS:'+str("% .2f" % orig_max_fos))
	sub3.text(0,z[7],'')
	sub3.text(0,z[8],'Optimized Design-', fontweight='bold')
	sub3.text(0,z[9],'Mass:'+str("% .3f" % opti_mass)+' kg')
	sub3.text(0,z[10],'Max. Disp.:'+str("% .2f" % opti_max_disp)+' mm')
	sub3.text(0,z[11],'Min. FoS:'+str("% .2f" % opti_min_fos))
	sub3.text(0,z[12],'Max. FoS:'+str("% .2f" % opti_max_fos))
	sub3.text(0,z[13],'')
	sub3.text(0,z[14],'Optimization Performance-',fontweight='bold')
	sub3.text(0,z[15],'Fitness : '+str("% .8f" % parent_fitness))
	sub3.text(0,z[16],'Time: '+str("% .2f" % elapsed_time)+" sec")
	sub3.text(0,z[17],'Iterations: '+str(i))
	sub3.text(0,z[18],'Successful Mutations: '+str(num_parents))
	
#main program
start_time = time.time()
gene_range = 50 / 1000 # mm converted to meters
max_mutations = 4 #per optimization cycle
fg_size = 8 #setup plot
fig = plt.figure(figsize=(fg_size*2,fg_size))
info_fig = plt.figure(figsize=(fg_size*2,fg_size))
#str_values = designLibrary_test.simple_truss1()
str_values = designLibrary_test.bridge_truss_1()
str_values_ini = copy.deepcopy(str_values)
variable_nodes, genes = gather_node_info()
parent = trussAnalysis_2D.run_analysis(str_values)
trussAnalysis_2D.plot_results(parent, fig,1, 2, 1, "Original Design")
dref = max(parent[10])
mref = parent[13]
parent_fitness, dfit, sfit, mfit = calculate_fitness(parent)
orig_min_fos = min(parent[14])
orig_max_fos = max(parent[14])
orig_max_disp = max(parent[10]) * 1000 #convert m to mm
orig_mass = parent[13]
i=0 ; x=[] ; y_parent=[] ; y_child=[] ; y_mass = [] ; y_disp = []
y_fos_min=[] ; y_fos_max=[] ; num_parents = 0
y_dfit=[] ; y_sfit=[] ; y_mfit=[]
while True:
	i+=1
	if i % 1000 == 0: # print status
		print(i, parent_fitness, child_fitness)
	genes, changes = mutate_genes()
	child = insert_genes()
	child = trussAnalysis_2D.run_analysis(child)
	child_fitness, dfit, sfit, mfit = calculate_fitness(child)
	x.append(i)
	y_parent.append(parent_fitness)
	y_child.append(child_fitness)
	y_mass.append(child[13])
	y_disp.append(max(child[10]) * 1000) #convert m to mm
	y_fos_min.append(min(child[14]))
	y_fos_max.append(max(child[14]))
	y_dfit.append(dfit)
	y_sfit.append(sfit)
	y_mfit.append(mfit)
	if child_fitness < parent_fitness:
		undo_mutate_genes()
	if child_fitness >= parent_fitness:
		parent_fitness = child_fitness
		parent = child[:]
		num_parents += 1
	if i == 20000 or keyboard.is_pressed('ctrl+q'):
#	if parent_fitness > 1.7 or keyboard.is_pressed('ctrl+q'):
		print('User halted the optimization.')
		print('solution fitness: ' + str("% .6f" % parent_fitness))
		print('iterations: '+ str(i))
		break
end_time = time.time()
elapsed_time = end_time - start_time
opti_min_fos = min(parent[14])
opti_max_fos = max(parent[14])
opti_max_disp = max(parent[10]) * 1000 #convert m to mm
opti_mass = parent[13]
print_results()
trussAnalysis_2D.plot_results(parent, fig, 1, 2, 2, "Optimized Design")
info_plotter()
trussAnalysis_2D.plot_results(parent, info_fig, 2, 3, 6, "Optimized Design")
fig.tight_layout()
plt.show()

#genes, changes = mutate_genes()
#child = insert_genes()
#undo_mutate_genes()

