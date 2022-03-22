#design library
import STL_reader

def simple_truss1():
	x = 1
	z = 1
	L1 = -1000
	A = 1 / 10000 # convert sq.cm to sq.m
	mat = 's' # input material code
	# NODE LISTS
	#	    coord  |   supports     | loads   |
	#        (x,y,z) |(0=fixed,1=free)| (x,y,z) |
	n1  = [ [ 0, 0,  0], [0, 0, 0], [0, 0, 0] ]
	n2  = [ [ x, 0,  z], [1, 0, 1], [0, 0, 0] ]
	n3  = [ [2*x, 0, 0], [1, 0, 1], [0, 0, L1] ]
	n4  = [ [3*x, 0, z], [1, 0, 1], [0, 0, 0] ]
	n5  = [ [4*x, 0, 0], [1, 0, 0], [0, 0, 0] ]
	# ELEMENT LISTS
	#    | node  | xcArea |  material |
	#    | pair  |(sq.cm) | (s, a, t) |
	E1 = [ [n1,n2], A, mat ]
	E2 = [ [n1,n3], A, mat ]
	E3 = [ [n2,n3], A, mat ]
	E4 = [ [n2,n4], A, mat ]
	E5 = [ [n3,n4], A, mat ]
	E6 = [ [n3,n5], A, mat ]
	E7 = [ [n4,n5], A, mat ]
	str_values = [ [E1,E2,E3,E4,E5,E6,E7] ]
	return str_values

#2D bridge structure, from COOK 4th ED., CH2, problem C2.5
def bridge_truss_1():
	A = 4 / 10000 # convert sq.cm to sq.m
	mat = 's' # input material code
	q = 1000  # N, applied loads, uniform on top surface
	xa = 0.0
	xb = 1.6
	xc = 3.2
	xd = 4.8
	xe = 6.4
	xf = 8.0
	xg = 9.6
	
	ya = 0.0
	yb = 0.2
	yc = 0.7
	yd = 1.7
	#coordinate (x,y) / supports (0=fixed, 1=free) / loads (x,y)
	n1  = [ [xa, 0, yd], [0, 0, 0], [0, 0, 0] ]  
	n2  = [ [xb, 0, yc], [1, 0, 1], [0, 0, 0] ]
	n3  = [ [xb, 0, yd], [1, 0, 1], [0, 0, -q] ]
	n4  = [ [xc, 0, yd], [1, 0, 1], [0, 0, -q] ]
	n5  = [ [xc, 0, yb], [1, 0, 1], [0, 0, 0] ]
	n6  = [ [xd, 0, yd], [1, 0, 1], [0, 0, -q] ]
	n7  = [ [xd, 0, ya], [1, 0, 1], [0, 0, 0] ]
	n8  = [ [xe, 0, yd], [1, 0, 1], [0, 0, -q] ]
	n9  = [ [xe, 0, yb], [1, 0, 1], [0, 0, 0] ]
	n10 = [ [xf, 0, yd], [1, 0, 1], [0, 0, -q] ]
	n11 = [ [xf, 0, yc], [1, 0, 1], [0, 0, 0] ]
	n12 = [ [xg, 0, yd], [1, 0, 0], [0, 0, 0] ]
	#-LIST ELEMENTS- 
	E1 = [ [n1,n2], A, mat ]
	E2 = [ [n1,n3], A, mat ]
	E3 = [ [n2,n3], A, mat ]
	E4 = [ [n3,n4], A, mat ]
	E5 = [ [n3,n5], A, mat ]
	E6 = [ [n2,n5], A, mat ]
	E7 = [ [n4,n5], A, mat ]
	E8 = [ [n4,n6], A, mat ]
	E9 = [ [n4,n7], A, mat ]
	E10 = [ [n5,n7], A, mat ]
	E11 = [ [n6,n7], A, mat ]
	E12 = [ [n6,n8], A, mat ]
	E13 = [ [n7,n8], A, mat ]
	E14 = [ [n7,n9], A, mat ]
	E15 = [ [n8,n9], A, mat ]
	E16 = [ [n8,n10], A, mat ]
	E17 = [ [n9,n10], A, mat ]
	E18 = [ [n9,n11], A, mat ]
	E19 = [ [n10,n11], A, mat ]
	E20 = [ [n10,n12], A, mat ]
	E21 = [ [n11,n12], A, mat ]
	str_values =[ [E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,E14,E15,E16,E17,E18,E19,E20,E21] ]
	return str_values
	
def cantilever_truss3D():
	print('input 3D truss structure parameters')
	x = 2
	y = 1
	z = 2
	L1 = -650 * 9.80665 #convert kg to N (kg*m/sec**2)
	DIA = 3 #cm
	A = (3.14159*DIA**2)/4 / 10000 # convert sq.cm to sq.m
	mat = 's' # input material code
	# NODE LISTS
			 #coord  |   supports     | loads   |
	#        (x,y,z) |(0=fixed,1=free)| (x,y,z) |
	n1  = [ [0, 0, 0], [0, 0, 0], [0, 0, 0] ]
	n2  = [ [0, -y, z],[0, 0, 0], [0, 0, 0] ]
	n3  = [ [0, y, z], [0, 0, 0], [0, 0, 0] ]
	n4  = [ [x, 0, z], [1, 1, 1], [0, 0, L1] ]

	# ELEMENT LISTS
	#    |nodes| xcArea | # material|
	#    |pair |(sq.cm) | #(s, a, w)|
	E1 = [ [n1,n4], A, mat ]
	E2 = [ [n2,n4], A, mat ]
	E3 = [ [n3,n4], A, mat ]

	str_values = [ [E1,E2,E3] ]
	return str_values
	
def cantileverBig_truss3D():
	print('input 3D truss structure parameters')
	x = 20
	y = 10
	z = 20
	L1 = -6500 * 9.80665 #convert kg to N (kg*m/sec**2)
	DIA = 25 #cm
	A = (3.14159*DIA**2)/4 / 10000 # convert sq.cm to sq.m
	mat = 's' # input material code
	# NODE LISTS
			 #coord  |   supports     | loads   |
	#        (x,y,z) |(0=fixed,1=free)| (x,y,z) |
	n1  = [ [0, 0, 0], [0, 0, 0], [0, 0, 0] ]
	n2  = [ [0, -y, z],[0, 0, 0], [0, 0, 0] ]
	n3  = [ [0, y, z], [0, 0, 0], [0, 0, 0] ]
	n4  = [ [x, 0, z], [1, 1, 1], [0, 0, L1] ]

	# ELEMENT LISTS
	#    |nodes| xcArea | # material|
	#    |pair |(sq.cm) | #(s, a, w)|
	E1 = [ [n1,n4], A, mat ]
	E2 = [ [n2,n4], A, mat ]
	E3 = [ [n3,n4], A, mat ]

	str_values = [ [E1,E2,E3] ]
	return str_values
	
def cantilever_2():
	x = 20
	y = 10
	z = 20
	L1 = -6500 * 9.80665 #convert kg to N (kg*m/sec**2)
	DIA = 12 #cm
	A = (3.14159*DIA**2)/4 / 10000 # convert sq.cm to sq.m
	mat = 's' # input material code
	# NODE LISTS
			 #coord  |   supports     | loads   |
	#        (x,y,z) |(0=fixed,1=free)| (x,y,z) |
	n1  = [ [0, 0, 0], [0, 0, 0], [0, 0, 0] ]
	n2  = [ [0, -y, z],[0, 0, 0], [0, 0, 0] ]
	n3  = [ [0, y, z], [0, 0, 0], [0, 0, 0] ]

	n4  = [ [x*(1/3), 0, z*(1/3)], [1, 1, 1], [0, 0, 0] ]
	n5  = [ [x*(1/3), -y*(2/3), z],[1, 1, 1], [0, 0, 0] ]
	n6  = [ [x*(1/3), y*(2/3), z], [1, 1, 1], [0, 0, 0] ]

	n7  = [ [x*(2/3), 0, z*(2/3)], [1, 1, 1], [0, 0, 0] ]
	n8  = [ [x*(2/3), -y*(1/3), z],[1, 1, 1], [0, 0, 0] ]
	n9  = [ [x*(2/3), y*(1/3), z], [1, 1, 1], [0, 0, 0] ]

	n10  = [ [x, 0, z], [1, 1, 1], [0, 0, L1] ]

	# ELEMENT LISTS
	#    |nodes| xcArea | # material|
	#    |pair |(sq.cm) | #(s, a, w)|
	E1 = [ [n1,n4], A, mat ]
	E2 = [ [n2,n5], A, mat ]
	E3 = [ [n3,n6], A, mat ]

	E4 = [ [n4,n5], A, mat ]
	E5 = [ [n5,n6], A, mat ]
	E6 = [ [n6,n4], A, mat ]

	E7 = [ [n4,n7], A, mat ]
	E8 = [ [n5,n8], A, mat ]
	E9 = [ [n6,n9], A, mat ]

	E10 = [ [n7,n8], A, mat ]
	E11 = [ [n8,n9], A, mat ]
	E12 = [ [n9,n7], A, mat ]

	E13 = [ [n7,n10], A, mat ]
	E14 = [ [n8,n10], A, mat ]
	E15 = [ [n9,n10], A, mat ]

	E16 = [ [n2,n4], A, mat ]
	E17 = [ [n3,n4], A, mat ]

	E18 = [ [n5,n7], A, mat ]
	E19 = [ [n6,n7], A, mat ]

	E20 = [ [n5,n9], A, mat ]
	E21 = [ [n6,n8], A, mat ]
	E22 = [ [n2,n6], A, mat ]
	E23 = [ [n3,n5], A, mat ]

	str_values = [ [E1,E2,E3,E4,E5,E6,E7,E8,E9,E10,E11,E12,E13,
	E14,E15,E16,E17,E18,E19,E20,E21,E22,E23] ]
	return str_values
	
def three_element_truss():
	filepath='C:/Users/Cameron/Desktop/'
	filename='3E_truss.stl'
	elements = STL_reader.STLreader(filepath, filename)
	str_values = []
	L1 = -650 * 9.80665 #convert kg to N (kg*m/sec**2)
	DIA = 5 #cm
	A = (3.14159*DIA**2)/4 / 10000 # convert sq.cm to sq.m
	mat = 's' # input material code
	fixd_sppt_locn = [ [0,0,0], [0,-1,2], [0,1,2] ]
	load_locn = [ [2,0,2] ]
	locn_tol = 0.05 #meters
	#keep only wanted elements that DO NOT span between fixed supports
	wanted_elements = []
	for element in elements:
#		print(element) #temporary diagnostic
		n1 = element[0][0][0] ; n2 = element[0][1][0]
		n1x = n1[0] ; n1y = n1[1] ; n1z = n1[2]
		n2x = n2[0] ; n2y = n2[1] ; n2z = n2[2]
		#for each element,
		#calculate the distance of each node to all support locations
		d1_list = []
		d2_list = []
		for locn in fixd_sppt_locn:
			sx = locn[0] ; sy = locn[1] ; sz = locn[2]
			d1 = ((n1x-sx)**2+(n1y-sy)**2+(n1z-sz)**2)**0.5
			d2 = ((n2x-sx)**2+(n2y-sy)**2+(n2z-sz)**2)**0.5
			d1_list.append(d1) ; d2_list.append(d2)
		#if both element nodes are at the same location as any of the
		#predefined support locations, DO NOT add the element to
		#the list of wanted elements. Otherwise, add element to the list.
#		print('min distance, node1 to any support:'+str(min(d1_list)))
#		print('min distance, node2 to any support:'+str(min(d2_list)))
		if min(d1_list) < locn_tol and min(d2_list) < locn_tol:
#			print('Element NOT added. \n')
			continue
		wanted_elements.append(element)
#		print('Element added. \n')
	elements = wanted_elements
#	print('list of wanted elements')
	
		#modify element data with support and load information.
	for element in elements:
		n1 = element[0][0][0] ; n2 = element[0][1][0]
		n1x = n1[0] ; n1y = n1[1] ; n1z = n1[2]
		n2x = n2[0] ; n2y = n2[1] ; n2z = n2[2]

		#apply support conditions
		for locn in fixd_sppt_locn:
			sx = locn[0] ; sy = locn[1] ; sz = locn[2]
			dist1 = ((n1x-sx)**2+(n1y-sy)**2+(n1z-sz)**2)**0.5
			dist2 = ((n2x-sx)**2+(n2y-sy)**2+(n2z-sz)**2)**0.5
			if dist1 < locn_tol: #if n1 is a fixed support location...
				#apply the appropriate fixed support condition to n1
				element[0][0][1] = [0,0,0]

			if dist2 < locn_tol: #if n2 is a fixed support location...
				#apply the appropriate fixed support condition to n2
				element[0][1][1] = [0,0,0]
				
		#apply loads
		for locn in load_locn:
			lx = locn[0] ; ly = locn[1] ; lz = locn[2]
			dist1 = ((n1x-lx)**2+(n1y-ly)**2+(n1z-lz)**2)**0.5
			dist2 = ((n2x-lx)**2+(n2y-ly)**2+(n2z-lz)**2)**0.5
			if dist1 < locn_tol:#if n1 is a load location...
				#apply the load at n1
				element[0][0][2] = [0, 0, L1]
			if dist2 < locn_tol:#if n2 is a load location...
				#apply the load at n2
				element[0][1][2] = [0, 0, L1]				
				
		element.append(A)
		element.append(mat)
		str_values.append(element)
	for e in str_values:
		print(e)

	return [str_values] # brackets are important!
#str_values = three_element_truss() #temporary
	
def transmissionTower1():
	print('importing simpleTower1 truss design...')
	L1 = -2000 * 9.80665 #convert kg to N (kg*m/sec**2)
	DIA = 6 #cm
	A = (3.14159*DIA**2)/4 / 10000 # convert sq.cm to sq.m
	mat = 's' # input material code
	filepath='C:/Users/Cameron/Desktop/'
	filename="simpleTower1.stl"
	elements = 	STL_reader.STLreader(filepath, filename)
	E_list = []
	str_values = []
	fixd_sppt_locn = [ [-1,1,0], [-1,-1,0], [1,-1,0],
	                   [1,1,0] ]
	load_locn = [ [4,0,12], [-4,0,12], [0,0,13] ]
	locn_tol = 0.01 # meters
	#no changes below this line
	
	#keep only wanted elements that DO NOT span between fixed supports
	wanted_elements = []
	for element in elements:
#		print(element) #temporary diagnostic
		n1 = element[0][0][0] ; n2 = element[0][1][0]
		n1x = n1[0] ; n1y = n1[1] ; n1z = n1[2]
		n2x = n2[0] ; n2y = n2[1] ; n2z = n2[2]
		#for each element,
		#calculate the distance of each node to all support locations
		d1_list = []
		d2_list = []
		for locn in fixd_sppt_locn:
			sx = locn[0] ; sy = locn[1] ; sz = locn[2]
			d1 = ((n1x-sx)**2+(n1y-sy)**2+(n1z-sz)**2)**0.5
			d2 = ((n2x-sx)**2+(n2y-sy)**2+(n2z-sz)**2)**0.5
			d1_list.append(d1) ; d2_list.append(d2)
		#if both element nodes are at the same location as any of the
		#predefined support locations, DO NOT add the element to
		#the list of wanted elements. Otherwise, add element to the list.
#		print('min distance, node1 to any support:'+str(min(d1_list)))
#		print('min distance, node2 to any support:'+str(min(d2_list)))
		if min(d1_list) < locn_tol and min(d2_list) < locn_tol:
#			print('Element NOT added. \n')
			continue
		wanted_elements.append(element)
#		print('Element added. \n')
	elements = wanted_elements
#	print('list of wanted elements')
	
		#modify element data with support and load information.
	for element in elements:
		n1 = element[0][0][0] ; n2 = element[0][1][0]
		n1x = n1[0] ; n1y = n1[1] ; n1z = n1[2]
		n2x = n2[0] ; n2y = n2[1] ; n2z = n2[2]

		#apply support conditions
		for locn in fixd_sppt_locn:
			sx = locn[0] ; sy = locn[1] ; sz = locn[2]
			dist1 = ((n1x-sx)**2+(n1y-sy)**2+(n1z-sz)**2)**0.5
			dist2 = ((n2x-sx)**2+(n2y-sy)**2+(n2z-sz)**2)**0.5
			if dist1 < locn_tol: #if n1 is a fixed support location...
				#apply the appropriate fixed support condition to n1
				element[0][0][1] = [0,0,0]

			if dist2 < locn_tol: #if n2 is a fixed support location...
				#apply the appropriate fixed support condition to n2
				element[0][1][1] = [0,0,0]
				
		#apply loads
		for locn in load_locn:
			lx = locn[0] ; ly = locn[1] ; lz = locn[2]
			dist1 = ((n1x-lx)**2+(n1y-ly)**2+(n1z-lz)**2)**0.5
			dist2 = ((n2x-lx)**2+(n2y-ly)**2+(n2z-lz)**2)**0.5
			if dist1 < locn_tol:#if n1 is a load location...
				#apply the load at n1
				element[0][0][2] = [0, 0, L1]
			if dist2 < locn_tol:#if n2 is a load location...
				#apply the load at n2
				element[0][1][2] = [0, 0, L1]				
				
		element.append(A)
		element.append(mat)
		str_values.append(element)
#	for e in str_values:
#		print(e)

	return [str_values] # brackets are important!


def spaceFrameRoof():
	L1 = -100 * 9.80665 #convert kg to N (kg*m/sec**2) about 1psi
	DIA = 2.91 #cm 2.91cm/FoS=2.5  ,  2.25cm/FoS=1.0
	A = (3.14159*DIA**2)/4 / 10000 # convert sq.cm to sq.m
	mat = 's' # input material code
	filepath="C:/Users/Cameron/Desktop/"
	filename="roof.stl"
	elements = STL_reader.STLreader(filepath, filename)
	E_list = [] ; str_values = []

	#define all support locations
	fixd_sppt_locn=[[.5,.5,-.473],[4.5,.5,-.473],[4.5,9.5,-.473],
	[.5,9.5,-.473]]
	slid_sppt_locn = [  ]
	supports = fixd_sppt_locn + slid_sppt_locn

	# calculate all load locations
	load_locn = [] ; num_x_rows = 6 ; num_y_rows = 11 ; z_coord = 0
	spacing = 1 #meters
	for i in range(num_x_rows):
		x_coord = spacing*i
		for j in range(num_y_rows):
			y_coord = spacing*j
			load_locn.append([x_coord, y_coord, z_coord])
	locn_tol = 0.01 # meters
	#no changes below this line
	
	#keep only wanted elements that DO NOT span between fixed supports
	wanted_elements = []
	for element in elements:
#		print(element) #temporary diagnostic
		n1 = element[0][0][0] ; n2 = element[0][1][0]
		n1x = n1[0] ; n1y = n1[1] ; n1z = n1[2]
		n2x = n2[0] ; n2y = n2[1] ; n2z = n2[2]
		#for each element,
		#calculate the distance of each node to all support locations
		d1_list = []
		d2_list = []
		for locn in fixd_sppt_locn:
			sx = locn[0] ; sy = locn[1] ; sz = locn[2]
			d1 = ((n1x-sx)**2+(n1y-sy)**2+(n1z-sz)**2)**0.5
			d2 = ((n2x-sx)**2+(n2y-sy)**2+(n2z-sz)**2)**0.5
			d1_list.append(d1) ; d2_list.append(d2)
		#if both element nodes are at the same location as any of the
		#predefined support locations, DO NOT add the element to
		#the list of wanted elements. Otherwise, add element to the list.
#		print('min distance, node1 to any support:'+str(min(d1_list)))
#		print('min distance, node2 to any support:'+str(min(d2_list)))
		if min(d1_list) < locn_tol and min(d2_list) < locn_tol:
#			print('Element NOT added. \n')
			continue
		wanted_elements.append(element)
#		print('Element added. \n')
	elements = wanted_elements
#	print('list of wanted elements')
	
		#modify element data with support and load information.
	for element in elements:
		n1 = element[0][0][0] ; n2 = element[0][1][0]
		n1x = n1[0] ; n1y = n1[1] ; n1z = n1[2]
		n2x = n2[0] ; n2y = n2[1] ; n2z = n2[2]

		#apply support conditions
		for locn in fixd_sppt_locn:
			sx = locn[0] ; sy = locn[1] ; sz = locn[2]
			dist1 = ((n1x-sx)**2+(n1y-sy)**2+(n1z-sz)**2)**0.5
			dist2 = ((n2x-sx)**2+(n2y-sy)**2+(n2z-sz)**2)**0.5
			if dist1 < locn_tol: #if n1 is a fixed support location...
				#apply the appropriate fixed support condition to n1
				element[0][0][1] = [0,0,0]

			if dist2 < locn_tol: #if n2 is a fixed support location...
				#apply the appropriate fixed support condition to n2
				element[0][1][1] = [0,0,0]
				
		#apply loads
		for locn in load_locn:
			lx = locn[0] ; ly = locn[1] ; lz = locn[2]
			dist1 = ((n1x-lx)**2+(n1y-ly)**2+(n1z-lz)**2)**0.5
			dist2 = ((n2x-lx)**2+(n2y-ly)**2+(n2z-lz)**2)**0.5
			if dist1 < locn_tol:#if n1 is a load location...
				#apply the load at n1
				element[0][0][2] = [0, 0, L1]
			if dist2 < locn_tol:#if n2 is a load location...
				#apply the load at n2
				element[0][1][2] = [0, 0, L1]				
				
		element.append(A)
		element.append(mat)
		str_values.append(element)
#	for e in str_values:
#		print(e)

	return [str_values] # brackets are important!
