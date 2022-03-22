#Import Libraries
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pylab as pl
import designLibrary_test

#Functions
def test_plot(str_values):
	print('plot undeformed structure to verify original design')
	fig = plt.figure() #setup the plot
	fig.subplots_adjust(bottom=0,top=1) #fit 3Dplot in window better
	sub = fig.add_subplot(1,1,1,projection="3d")
	X = []  ;  Y = []  ;  Z = []
	#retrieve node/element into from data
	ele_values = str_values[0]
	for ele_prop in ele_values:
		x = []  ;  y = []  ;  z = []
		#retrieve x,z,y components of nodal coordinates
		#x: horizontal, y: out of plane, z: vertical
		n1x = ele_prop[0][0][0][0]
		n1y = ele_prop[0][0][0][1]
		n1z = ele_prop[0][0][0][2]
		n2x = ele_prop[0][1][0][0]
		n2y = ele_prop[0][1][0][1]
		n2z = ele_prop[0][1][0][2]
		x.append(n1x)  ;  X.append(n1x)
		x.append(n2x)  ;  X.append(n2x)
		y.append(n1y)  ;  Y.append(n1y)
		y.append(n2y)  ;  Y.append(n2y)
		z.append(n1z)  ;  Z.append(n1z)
		z.append(n2z)  ;  Z.append(n2z)
		sub.plot(x,y,z, color = 'blue', marker = 'o', markersize = 3)
	# Create cubic bounding box to simulate equal aspect ratio
	box_size = max( max(X)-min(X), max(Y)-min(Y), max(Z)-min(Z) )
	a = -( box_size/2 - (max(X) + min(X))/2 )
	b =    box_size/2 + (max(X) + min(X))/2
	c = -( box_size/2 - (max(Y) + min(Y))/2 )
	d =    box_size/2 + (max(Y) + min(Y))/2
	e = -( box_size/2 - (max(Z) + min(Z))/2 )
	f =    box_size/2 + (max(Z) + min(Z))/2
	x_b = [a,b]
	y_b = [d,c]
	z_b = [f,e]
	sub.plot(x_b, y_b, z_b, linewidth=0, color='red')
	#after all elements have been added to the plot, display the figure
	plt.show()
	
def get_material_properties(material):
  if material == 's':
    #structural steel, ASTM A36
    E = 200 * 10**9 #Pa
    s_max = 250 * 10**6 #Pa yield stress
    density = 7850 #Kg/cu.m
  if material == 'a':
    #standard 6061-T6 alloy
    E = 102 * 10**9 #Pa
    s_max = 276 * 10**6 #Pa yield stress
    density = 2700 #Kg/cu.m
  if material == 'w':
    #wood, douglas-fir
    E = 10.8 * 10**9 #Pa
    s_max = 21.6 * 10**6 #Pa yield stress
    density = 526 #Kg/cu.m
  if material == 't':
    #based on most widely used 6Al-4V alloy
    E = 102 * 10**9 #Pa
    s_max = 340 * 10**6 #Pa yield stress
    density = 4510 #Kg/cu.m
  return E, s_max, density
  
def calculate_element_values(str_values):
	#print('calculate element values')
	ele_values = str_values[0]
	for ele_props in ele_values:
		#calculate element length
		n1x = ele_props[0][0][0][0]
		n1y = ele_props[0][0][0][1]
		n1z = ele_props[0][0][0][2]
		n2x = ele_props[0][1][0][0]
		n2y = ele_props[0][1][0][1]
		n2z = ele_props[0][1][0][2]
		ele_len = ( (n2x-n1x)**2+(n2y-n1y)**2+(n2z-n1z)**2 )**(1/2)
		#get element material properties
		material = ele_props[2]
		E, s_max, density = get_material_properties(material)
		#convert xc_area value from sq.cm to sq.m
		ele_xc_area = ele_props[1]
		#calculate element mass
		ele_mass = ele_len * ele_xc_area * density #Kg
		#add values to element properties list
		ele_props.append(ele_len)
		ele_props.append(E)
		ele_props.append(s_max)
		ele_props.append(ele_mass)
	ele_values = str_values[0]
#	for ele_props in ele_values:
#		print(ele_props)
		
def add_self_weight_loads(str_values):
	#unless vertically fixed, each node gets 1/2 the weight of every
	#element that is attached to it. Note wt loads are negative (-z dir)
	#print("add self-weight loads")
	grav = 9.807 #m/sec^2
	ele_wt_total = 0
	added_load = 0
	lost_load = 0
	ele_values = str_values[0]
	for ele_props in ele_values:
		ele_mass = ele_props[6]
		ele_wt = ele_mass * grav #Newtons
		ele_wt_total += ele_wt
		#if the 1st element node is free (=1) in the vertical direction,
		#increment the vertical load on that node by 1/2 element wt.
		supt1_z = ele_props[0][0][1][2]	
		if supt1_z == 1:
			ele_props[0][0][2][2] += -ele_wt/2
			added_load += ele_wt/2
		else:
			lost_load += ele_wt/2
			
		#if the 2nd element node is free (=1) in the vertical direction,
		#increment the vertical load on that node by 1/2 element wt.
		supt2_z = ele_props[0][1][1][2]
		if supt2_z == 1:
			ele_props[0][1][2][2] += -ele_wt/2
			added_load += ele_wt/2
		else:
			lost_load += ele_wt/2
#	print("\ntotal structure wt: ", ele_wt_total)
#	print("added load: ", added_load)
#	print("lost load: ", lost_load)
#	print("added + lost load: ", added_load + lost_load)
#	for ele_props in ele_values:
#		print(ele_props)
		
def calculate_element_stiffness_matrices(str_values):
	#print('calculate element stiffness matrices')
	ele_values = str_values[0]
	for ele_props in ele_values:
		A = ele_props[1]
		E = ele_props[4]
		L = ele_props[3]
		n1x = ele_props[0][0][0][0]
		n1y = ele_props[0][0][0][1]
		n1z = ele_props[0][0][0][2]
		n2x = ele_props[0][1][0][0]
		n2y = ele_props[0][1][0][1]
		n2z = ele_props[0][1][0][2]
		Vx = n2x - n1x
		Vy = n2y - n1y
		Vz = n2z - n1z
		Lv = ( Vx**2 + Vy**2 + Vz**2 ) **(1/2)
		l = Vx / Lv
		m = Vy / Lv
		n = Vz / Lv
		T = np.array([[l**2,  l*m,  l*n,  -(l**2),  -l*m,  -l*n],
				  [m*l,  m**2,  m*n,  -m*l,  -(m**2),  -m*n],
				  [n*l,  n*m,  n**2,  -n*l,  -n*m,  -(n**2)],
				  [-(l**2),  -l*m,  -l*n,  l**2,  l*m,  l*n],
				  [-m*l,  -(m**2),  -m*n,  m*l,  m**2,  m*n],
				  [-n*l,  -n*m,  -(n**2),  n*l,  n*m,  n**2]])
		k = (E*A/L)
		esm = k * T
		ele_props.append(esm)
#	np.set_printoptions(precision = 2, linewidth=150)
#	ele_values = str_values[0]
#	for ele_props in ele_values:
#		print(ele_props[7])
		
def generate_global_stiffness_matrix(str_values):
	#print('generate global stiffness matrix')
#	print('make a list of all nodes')
	node_list = []
	ele_values = str_values[0]
	for ele_props in ele_values:
		n1x = ele_props[0][0][0][0]
		n1y = ele_props[0][0][0][1]		
		n1z = ele_props[0][0][0][2]
		n2x = ele_props[0][1][0][0]
		n2y = ele_props[0][1][0][1]		
		n2z = ele_props[0][1][0][2]
		node1 = [n1x, n1y, n1z]
		node2 = [n2x, n2y, n2z]
		if node1 not in node_list:
			node_list.append(node1)
		if node2 not in node_list:
			node_list.append(node2)
	str_values.append(node_list)
#	print('make a list of global DOFs')
	dof_list = []
	for i in range(0, len(node_list)*3, 3):
		dof_3D = [i,i+1,i+2]
		dof_list.append(dof_3D)
	str_values.append(dof_list)
#	print(dof_list)
#	print('assign DOFs to elements')
	global_dof_list = []
	for ele_props in ele_values:
		n1x = ele_props[0][0][0][0]
		n1y = ele_props[0][0][0][1]
		n1z = ele_props[0][0][0][2]
		n2x = ele_props[0][1][0][0]
		n2y = ele_props[0][1][0][1]
		n2z = ele_props[0][1][0][2]
		node1_3D = [n1x, n1y, n1z]
		node2_3D = [n2x, n2y, n2z]
		for i in range(len(node_list)):
			if node_list[i] == node1_3D:
				dof1 = dof_list[i]
			if node_list[i] == node2_3D:
				dof2 = dof_list[i]
		dof_pair = [dof1, dof2]
		global_dof_list.append(dof_pair)
	str_values.append(global_dof_list)
#	print(global_dof_list)
#	print('destination DOF')
	destination_dof = []
	for pair in global_dof_list:
#		print(pair)
		dof_group = []
		for item in pair:
#			print(item)
			for number in item:
				dof_group.append(number)
		destination_dof.append(dof_group)
	str_values.append(destination_dof)
#	print(destination_dof)
	#make a zeroes matrix to match the size of the gsm
	#each node has 3 degrees of freedom...
	size = ( 3*len(node_list), 3*len(node_list) )
	gsm = np.zeros(size)
#	print(gsm)
	#map ESMs to the GSM !!!
	for i in range(0, len(ele_values)):
		esm = ele_values[i][7]
		dd = destination_dof[i]
		i = 0
		for x in dd:
			j = 0
			for y in dd:
				#i,j is the ESM value
				#x,y is the GSM destination
				gsm[x][y] += esm[i][j]
				j += 1
			i += 1
#	print('apply boundary conditions')
#make a list of support conditions for each global DOF
	supports = []
	support_condition = []
	for i in range(0, len(ele_values)):
		s_n1x = ele_values[i][0][0][1][0]
		s_n1y = ele_values[i][0][0][1][1]
		s_n1z = ele_values[i][0][0][1][2]
		s_n2x = ele_values[i][0][1][1][0]
		s_n2y = ele_values[i][0][1][1][1]
		s_n2z = ele_values[i][0][1][1][2]
		support_condition.append([s_n1x, s_n1y, s_n1z, s_n2x, s_n2y, s_n2z])
	used_dof = []
	for dof in range(0, len(node_list*3)):
		for j in range(0, len(ele_values)):
			for k in range(0,len(destination_dof[0])):
				if destination_dof[j][k] == dof and dof not in used_dof:
					supports.append(support_condition[j][k])
					used_dof.append(dof)
#	print(supports)
	#modify GSM with the boundary conditions...
	for i in range(0, len(supports)):
		support = supports[i]
		if support == 0: # if DOF is fixed...
			gsm[i] = 0   # zero out row
			gsm[:,i] = 0 # zero out column
			gsm[i,i] = 1 # make the diagonal position equal 1
	str_values.append(gsm)
#	print('global stiffness matrix with boundary conditions applied')
#	print(gsm)
#		str_values.append(gsm)
#	print(gsm)

def create_loads_vector(str_values):
	ele_values = str_values[0]
	node_list = str_values[1]
	destination_code = str_values[4]
	#create the loads vector
	#print('create the loads vector')
	loads = []
	load_condition = []
	for i in range(0, len(ele_values)):
		L_n1x = ele_values[i][0][0][2][0]
		L_n1y = ele_values[i][0][0][2][1]
		L_n1z = ele_values[i][0][0][2][2]
		L_n2x = ele_values[i][0][1][2][0]
		L_n2y = ele_values[i][0][1][2][1]
		L_n2z = ele_values[i][0][1][2][2]
		load_condition.append([L_n1x,L_n1y,L_n1z,L_n2x,L_n2y,L_n2z])
	used_i = []
	for i in range(0, len(node_list*3)):
		for j in range(0, len(ele_values)):
			for k in range(0,len(destination_code[0])):
				if destination_code[j][k] == i and i not in used_i:
					loads.append(load_condition[j][k])
					used_i.append(i)
	str_values.append(loads)
#	print(loads)
	
def calculate_displacements(str_values):
	#print('calculate displacements')
	gsm = str_values[5]
	loads = str_values[6]
	nodal_disp = np.linalg.solve(gsm,loads)
	str_values.append(nodal_disp)
#	print(nodal_disp)
	
def calculate_deformed_element_length(str_values):
	#print('calculate new nodal positions') #change to element length?
	node_list = str_values[1]
	nodal_disp = str_values[7]
	displaced_nodes = []
	for i in range(0, len(node_list)):
		node_x = node_list[i][0]
		node_y = node_list[i][1]
		node_z = node_list[i][2]
		disp_x = nodal_disp[3*i]
		disp_y = nodal_disp[3*i+1]
		disp_z = nodal_disp[3*i+2]
		disp_node = [(node_x + disp_x),(node_y + disp_y),(node_z + disp_z)]
		displaced_nodes.append(disp_node)
	str_values.append(displaced_nodes)
	ele_values = str_values[0]
	ele_len_deformed_list = []
	for ele_props in ele_values:
		n1x = ele_props[0][0][0][0]
		n1y = ele_props[0][0][0][1]
		n1z = ele_props[0][0][0][2]
		n2x = ele_props[0][1][0][0]
		n2y = ele_props[0][1][0][1]
		n2z = ele_props[0][1][0][2]
		node1_3D = [n1x, n1y, n1z]
		node2_3D = [n2x, n2y, n2z]
		for i in range(0, len(node_list)):
			if node_list[i] == node1_3D:
				n1d = displaced_nodes[i]
			if node_list[i] == node2_3D:
				n2d = displaced_nodes[i]
		n1dx = n1d[0]  ;  n1dy = n1d[1] ;  n1dz = n1d[2]
		n2dx = n2d[0]  ;  n2dy = n2d[1] ;  n2dz = n2d[2]
		ele_len_deformed = ( (n2dx-n1dx)**2+(n2dy-n1dy)**2+(n2dz-n1dz)**2 )**(1/2)
		ele_len_deformed_list.append(ele_len_deformed)
	str_values.append(ele_len_deformed_list)

	# calculate magnitude change of nodal position
	node_position_magnitude_change = []
	disp_nodes = str_values[8]
	for i in range(0, len(node_list)):
		orig_pos = node_list[i]
		disp_pos = disp_nodes[i]
		n1dx = orig_pos[0]; n1dy = orig_pos[1]; n1dz = orig_pos[2]
		n2dx = disp_pos[0]; n2dy = disp_pos[1]; n2dz = disp_pos[2]
		move_dist = ( (n2dx-n1dx)**2+(n2dy-n1dy)**2+(n2dz-n1dz)**2 )**(1/2)
		node_position_magnitude_change.append(move_dist)
	str_values.append(node_position_magnitude_change)
#	print(node_position_magnitude_change)

def calculate_strain_and_stress(str_values):
	#print('calculate element  strains & stresses')
	ele_strain = []
	ele_stress = []
	ele_len_deformed_list = str_values[9]
	ele_values = str_values[0]
	for i in range(0, len(ele_values)):
		ele_len_original = ele_values[i][3]
		modulus = ele_values[i][4]
		ele_len_deformed = ele_len_deformed_list[i]
		strain = (ele_len_deformed - ele_len_original)/ele_len_original
		stress = strain * modulus
		ele_strain.append(strain)
		ele_stress.append(stress)
	str_values.append(ele_stress)
	str_values.append(ele_strain)
#	print(ele_strain)
#	print(ele_stress)
	
def calculate_other_results(str_values):
	#print('calculate other results')
	#total structure mass
	total_mass = 0
	ele_values = str_values[0]
	for ele_props in ele_values:
		ele_mass = ele_props[6]
		total_mass += ele_mass
	str_values.append(total_mass)
#	print(total_mass)
	#element safety factor, tension failure or buckling failure
	failure_mode = []
	ele_factor_of_safety = []
	ele_stress = str_values[11]
	for i in range(len(ele_values)):
		ele_calcd_stress = ele_stress[i]
		ele_mod = ele_values[i][4] #modulus of elasticity
		ele_max_allow_stress = ele_values[i][5] #tension
		if ele_calcd_stress >= 0: # element in tension
			fos_t = ele_max_allow_stress / ele_calcd_stress
			ele_factor_of_safety.append(fos_t)
			failure_mode.append('t')

		if ele_calcd_stress < 0: #element in compression
			K = 1 # column effective length factor
			 # = 1 for pinned-pinned, =0.5 for fixed-fixed
			pi = 3.14159
			ele_xc_area = ele_values[i][1]
			ele_r = (ele_xc_area/pi)**.5
			I_x_circle = (pi/4)*ele_r**4
			l_e= ele_values[i][3]*K #effective element length
			r_gyr = (I_x_circle/ele_xc_area)**0.5 #radius of gyration
			#calculate element critical buckling stress
			ele_cb_stress = -(pi**2*ele_mod)/((l_e/r_gyr)**2)
			#calculate compressive yield stress
			ele_cy_stress = -ele_max_allow_stress / ele_calcd_stress
			#use the smaller of the two critical stresses to calc
			#factor of safety for compressively loaded elements
			if ele_cb_stress <= ele_cy_stress:
				fos = ele_cb_stress / ele_calcd_stress
				failure_mode.append('b')
			else:
				fos = ele_cy_stress / ele_calcd_stress
				failure_mode.append('c')
			ele_factor_of_safety.append(fos)
		str_values.append(ele_factor_of_safety)
		str_values.append(failure_mode)
		
	#calculate strain energy of each element - old version
#	strain_energy = []
#	ele_area = []
#	ele_len = []
#	ele_E = []
#	for ele_props in ele_values:
#		e_area = ele_props[1] ; ele_area.append(e_area)
#		e_len = ele_props[3] ; ele_len.append(e_len)
#		e_E = ele_props[4] ; ele_E.append(e_E)
#	for i in range(0, len(ele_values)):
#		s_nrg = 0.5*ele_area[i]*ele_len[i]*(ele_stress[i]**2)/ ele_E[i]
#		strain_energy.append(s_nrg)
#	str_values.append(strain_energy) #str_values[16]
	
	#calculate strain energy of each element
	strain_energy = []
	ele_stress = str_values[11]
	i=0
	for ele_props in ele_values:
		e_A = ele_props[1] #element cross sectional area
		e_L = ele_props[3] #element length
		e_E = ele_props[4] #element elastic modulus
		s_nrg = 0.5*e_A*e_L*(ele_stress[i]**2) / e_E #strain energy
		strain_energy.append(s_nrg)
		i+=1
		str_values.append(strain_energy) #str_values[16]

#	print(str_values[13])
#	print(str_values[14])
#	print(str_values[15])

def plot_results(str_values, fig, rows, cols, subplot_num, title):
	#print('plot results of analysis')
	#setup the figure
	#fig = plt.figure(figsize=(9,9)) #setup the plot
	sub = fig.add_subplot(rows,cols,subplot_num,projection="3d")
	sub.view_init(elev=25, azim=280)
	sub.dist = 9
	sub.set_title(title)
	sub.set_xlabel('X', fontweight='bold')
	sub.set_ylabel('Y', fontweight='bold')
	sub.set_zlabel('Z', fontweight='bold')
	sub.set_facecolor('lightgrey')
	#retrieve node/element into from data
	ele_values = str_values[0]
	node_list = str_values[1]
	nodal_disp = str_values[7]
	displaced_nodes = str_values[8]
	ele_stress = str_values[11]
	total_mass = str_values[13]
	factor_of_safety = str_values[14]
	failure_mode = str_values[15]
	# plot original structure
	X = []  ;  Y = []  ;  Z = []
	for ele_prop in ele_values:
		x = []  ;  y = []  ;  z = []
		# retrieve x,y,z components of nodal coordinates
		# x: horizontal, y: out of plane, z: vertical
		n1x = ele_prop[0][0][0][0]
		n1y = ele_prop[0][0][0][1]
		n1z = ele_prop[0][0][0][2]
		n2y = ele_prop[0][1][0][1]
		n2x = ele_prop[0][1][0][0]
		n2z = ele_prop[0][1][0][2]
		x.append(n1x)  ;  X.append(n1x)
		x.append(n2x)  ;  X.append(n2x)
		y.append(n1y)  ;  Y.append(n1y)
		y.append(n2y)  ;  Y.append(n2y)
		z.append(n1z)  ;  Z.append(n1z)
		z.append(n2z)  ;  Z.append(n2z)
		sub.plot(x,y,z, color = 'black', marker = 'o',
		markersize = 1, linewidth=.5, linestyle='dashed')
	# Create cubic bounding box to simulate equal aspect ratio
	box_size = max( max(X)-min(X), max(Y)-min(Y), max(Z)-min(Z) )
	a = -( box_size/2 - (max(X) + min(X))/2 )
	b =    box_size/2 + (max(X) + min(X))/2
	c = -( box_size/2 - (max(Y) + min(Y))/2 )
	d =    box_size/2 + (max(Y) + min(Y))/2
	e = -( box_size/2 - (max(Z) + min(Z))/2 )
	f =    box_size/2 + (max(Z) + min(Z))/2
	x_b = [a,b]
	y_b = [d,c]
	z_b = [f,e]
	sub.plot(x_b, y_b, z_b, linewidth=0)
	#plot displaced structure
	#set up the colormap!
	n = 40 #number of color levels in the color linspace
	c_factor = 1.15 #max/min multiplier so colors aren’t so dark
	tension_colors = pl.cm.Reds(np.linspace(0,1,n))
	compression_colors = pl.cm.Blues(np.linspace(0,1,n))
	# create a list of displaced nodes- not scaled and scaled
	max_disp = (max(str_values[10]) * 1000)	#convert m to mm
	disp_scale = 25/max_disp
	displaced_nodes_scaled = []
	for i in range(0, len(node_list)):
		node_x = node_list[i][0]
		node_y = node_list[i][1]
		node_z = node_list[i][2]
		disp_x_s = nodal_disp[3*i] * disp_scale #mod, changed 2 to 3
		disp_y_s = nodal_disp[3*i+1] * disp_scale
		disp_z_s = nodal_disp[3*i+2] * disp_scale
		disp_node_s = [(node_x + disp_x_s),(node_y + disp_y_s),(node_z + disp_z_s)]
		displaced_nodes_scaled.append(disp_node_s)
		
	s_index = 0
	s_max = max(ele_stress)
	s_min = min(ele_stress)
	for ele_props in ele_values:
		stress = ele_stress[s_index]
		fm = failure_mode[s_index]
		fos = factor_of_safety[s_index]
		s_index += 1
		#calculate color_num
		if stress > 0:
			color_num = int(round(n*((stress)/(s_max*c_factor))))
		if stress < 0:
			color_num = int(round(n*((stress)/(s_min*c_factor))))
		#condition color_num
		if color_num >= n: # high limit
			color_num = int(n-1)
		if color_num < 0: # low limit
			color_num = int(1)
		#configure element color, tension in red, compression in blue.
		if stress > 0:
			select_color = tension_colors[color_num]
		if stress < 0:
			select_color = compression_colors[color_num]
		n1x = ele_props[0][0][0][0]
		n1y = ele_props[0][0][0][1]
		n1z = ele_props[0][0][0][2]
		n2x = ele_props[0][1][0][0]
		n2y = ele_props[0][1][0][1]
		n2z = ele_props[0][1][0][2]
		node1_3D = [n1x, n1y, n1z]
		node2_3D = [n2x, n2y, n2z]
		for i in range(0, len(node_list)):
			x_d = []  ;  y_d = []  ;  z_d = []
			if node_list[i] == node1_3D:
				n1_d = displaced_nodes_scaled[i]
			if node_list[i] == node2_3D:
				n2_d = displaced_nodes_scaled[i]
		n1_dx = n1_d[0]  ;  n2_dx = n2_d[0]
		n1_dy = n1_d[1]  ;  n2_dy = n2_d[1]
		n1_dz = n1_d[2]  ;  n2_dz = n2_d[2]
		x_d.append(n1_dx) 
		x_d.append(n2_dx)
		y_d.append(n1_dy) 
		y_d.append(n2_dy) 			
		z_d.append(n1_dz) 
		z_d.append(n2_dz)
		#calculate location for element label
		x_ave = sum(x_d)/len(x_d)
		y_ave = sum(y_d)/len(y_d)
		z_ave = sum(z_d)/len(z_d)
		sub.plot(x_d,y_d,z_d, color = select_color,
		linewidth=3, linestyle='solid')
		#place a label on each element
		s_mpa = stress / 10**6
#		sub.text(x_ave, y_ave, z_ave,
#		str( "%.1f" % s_mpa +'MPa|FS:'+"%.2f" % fos+'('+fm+')'), fontsize=8.5,
#		fontweight='bold', color='black')
		#draw a dot at each node location
		sub.scatter(x_d,y_d,z_d, s=20, c='black', marker='o')
	#calculate & plot magnitude of change of nodal position
#	for i in range(0, len(node_list)):
#		orig_pos = node_list[i]
#		disp_pos = displaced_nodes[i]
#		n1dx = orig_pos[0]; n1dy = orig_pos[1]; n1dz = orig_pos[2]
#		n2dx = disp_pos[0]; n2dy = disp_pos[1]; n2dz = disp_pos[2]
#		move_dist = ( (n2dx-n1dx)**2+(n2dy-n1dy)**2+(n2dz-n1dz)**2 )**(1/2)
#		move_dist_mm = move_dist * 1000 # convert m to mm
#		sub.text(n2dx, n2dy, n2dz, str( "%.2f" % move_dist_mm + " mm"),
#		fontsize=8, fontweight='bold', color='darkblue')
	#display total structure mass
	sub.text(1,0,1.5,'Total Mass:'+str("% .3f" % total_mass)+' kg')
	#plt.show() #disable for multiplot, enable for single analysis
	
	
#Main Program
def run_analysis(str_values):
	#test_plot(str_values)
	calculate_element_values(str_values)
	add_self_weight_loads(str_values)
	calculate_element_stiffness_matrices(str_values)
	generate_global_stiffness_matrix(str_values)
	create_loads_vector(str_values)
	calculate_displacements(str_values)
	calculate_deformed_element_length(str_values)
	calculate_strain_and_stress(str_values)
	calculate_other_results(str_values)
	#plot_results(str_values, fig, 1, 1, 1, "3D Truss Analysis")
	return(str_values)


#fig = plt.figure(figsize=(8,8)) #temp
#str_values = designLibrary_test.three_element_truss() #temp
#str_values = designLibrary_test.transmissionTower1() #temp
#run_analysis(str_values) #temp

#print("strain energy (U) of each element")
#print(str_values[16])
#print(len(str_values[16]))

'''
#print ESMs
ele_props = str_values[0]
for prop in ele_props:
	print(np.array2string(prop[7], max_line_width=np.inf, precision = 0))
	print('\n')

#print GSM
print('node_list') ; print(str_values[1])
print('dof_list') ; print(str_values[2])
print('global_dof_list') ; print(str_values[3])
print('destination_dof') ; print(str_values[4])
print('gsm')
print(np.array2string(str_values[5], max_line_width=np.inf, precision = 3))

#print loads vector
print('loads') ; print(str_values[6])

#print displacements
print('displacements')
print(np.array2string(str_values[7], max_line_width=np.inf, precision = 6))

#print disp nodes, element lengths, node pos mag change
print('displaced nodes') ; print(str_values[8])
print('deformed element length') ; print(str_values[9])
print('node position magnitude change') ; print(str_values[10])
'''
