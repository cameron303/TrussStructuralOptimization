#define filepath and filename
#filepath='C:/Users/Cameron/Desktop/'
#filename='3E_truss.stl'
def STLreader(filepath, filename):
	print('Extracting information from the STL file...')
	f = open(filepath + filename, 'r')
	data = f.read()
	data = data.split()
	all_node_list = []
	unique_node_list = []
	recorded_nodes = []
	#make lists of 1.)all nodes, and 2.)unique nodes
	for i in range(0, len(data)):
		word = data[i]
		if word == 'vertex':
			v_x = float(data[i+1])
			v_y = float(data[i+2])
			v_z = float(data[i+3])
			node = [v_x, v_y, v_z]
			all_node_list.append(node)
			if node not in recorded_nodes:
				recorded_nodes.append(node)
				unique_node_list.append([node, [1,1,1],[0,0,0]])			

	# make a list of unique elements 
	element_list = []
	for i in range(0, int(len(all_node_list)),3):
		va = all_node_list[i]
		vb = all_node_list[i+1]
		vc = all_node_list[i+2]
		if [va,vb] not in element_list and [vb,va] not in element_list:
			element_list.append([va,vb])
		if [vb,vc] not in element_list and [vc,vb] not in element_list:
			element_list.append([vb,vc])
		if [vc,va] not in element_list and [va,vc] not in element_list:
			element_list.append([vc,va])
			
	#make a new formatted, referenced list of elements
	new_element_list = []
	for element in element_list:
		new_element = []
		node1 = element[0]
		node2 = element[1]
		node_pair = []
		for node in unique_node_list:
			if node[0] == node1:
				node_pair.append(node)
			if node[0] == node2:
				node_pair.append(node)
		new_element.append(node_pair)
		new_element_list.append(new_element)
	f.close()
	return new_element_list
'''
	#temporary diagnostic print statement
	print("\n all_node_list")
	for n in all_node_list:
		print(n)
	print("\n unique_node_list")
	for n in unique_node_list:
		print(n)
	print("\n unique_element_list")
	for e in element_list:
		print(e)
	print("\n formatted_element_list")		
	for e in new_element_list:
		print(e)


STLreader(filepath, filename) #temporary function call
'''
