'''
@ The script is developed by Qingbin Liu under the supervision of Jie Liu
E-mail to Qingbin Liu: liuqingb@mail2.sysu.edu.cn; liuqingb@pku.edu.cn
       to discuss technical details.
--------------------------------------------------------------------------------
The script is to convert the data from Abaqus output database format 
to vtk file format for parallel visualization
Python script performs the following three major steps: 
	1) reading Abaqus output file according to the architecture of ODB; 
	2) data decomposition for parallel visualization; 
	3) writing VTK for Paraview.
'''

#import necessary modules to handle Abaqus output database, files and string
from odbAccess import *
from textRepr import *
from string import *
from time import *
import numpy as np
import os
def ConvertOdb2VtkP(Odbpath = ' ',Odbname = ' ',vtkpath = ' ', BeginFrame = ' ', EndFrame =' ', Steps =' ', Instances = ' '):
	
	starttime = time()
	
	#get odb file's path
	odb_path = Odbpath
	#get odb file name
	odbname = Odbname
	#get the output files' path
	vtk_path = vtkpath

	#get the frame
	#input_frame = range(int(BeginFrame),int(EndFrame)+1)
	input_frame = range(0,2000,100)
	#get the step
	input_step = Steps.split(",")
	#get the instance
	input_instance = Instances.split(",")

	#display the read parameter
	print "Basic Information:"
	print "Model:",odbname
	print "Convert frames: ",input_frame[0]," to ",input_frame[-1]
	print "Step & Instance : ",str(input_step),", ",str(input_instance)
	
	#open an ODB ( Abaqus output database )
	odb = openOdb(os.path.join(odb_path,odbname)+'.odb',readOnly=True)
	print "ODB opened"
	
	#access geometry and topology information ( odb->rootAssembly->instances->(nodes, elements) )
	rootassembly = odb.rootAssembly
	instance = rootassembly.instances
	#access attribute information
	step = odb.steps
	#get instance & step information : Quantity and all names
	allinstancestr = str(instance)
	autoins = allinstancestr.split("'")
	inslen = len(autoins)/4
	instance_N = range(0,inslen)
	allstepstr = str(step)
	autostep = allstepstr.split("'")
	steplen = len(autostep)/4
	step_N = range(0,steplen)
	
	for i in input_step:
	    if(steplen < int(i)):
			print "input step exceeds the range of steps"
			os._exit(0)
	for i in input_instance:
	    if(inslen < int(i)):
			print "input instance exceeds the range of instances"
			os._exit(0)
		

	#step cycle
	for step_i in input_step:
		n = int(step_i)*4+1
		stepname = autostep[n]
		print "Step: ",stepname
		#instance cycle
		for ins_i in input_instance:
			n = int(ins_i)*4+1
			instancename = autoins[n]
			print "Instance: ",instancename
			
			#access nodes & elements
			node = instance[instancename].nodes
			element = instance[instancename].elements
			n_nodes = len(node)
			#access attribute(fieldOutputs) information
			frame = step[stepname].frames
			
					
			
			#match nodes' label and its order in sequence (for empty nodes in tetra mesh)
			MLN = node[n_nodes-1].label
			TOTAL=[]
			#read node in sequence, and get the largest label of node(non-empty) 
			#MLN is the max label of nodeset
			for i in node:
				TOTAL.append(i.label)
				if(i.label > MLN):
					MLN = i.label
			#match (the key)
			L=[]
			n = 0
			for i in range(MLN): 
				L.append(0)
			for i in TOTAL:
				L[i-1] = n
				n += 1

			#match element in sequence, and get the largest label of element
			MLE=element[0].label
						#read node in sequence, and get the largest label of node(non-empty) 
			#MLN is the max label of nodeset
			for i in element:
				if(i.label > MLE):
					MLE = i.label
			#match (the key)
			L_E=[]
			n = 0
			n_elements=0
			TOTAL_E=[]
			n_elenode=[]
			L_E=np.zeros(MLE)-1
			element_mark=np.zeros(MLE)
			for i in element:
			    L_E[i.label-1] = n
			    n += 1
			    if (i.type=='CPE3'):
			        TOTAL_E.append(i.label)
			        n_elenode.append(3)
			        n_elements+=1
			        element_mark[i.label-1]=1
			    elif (i.type=='CPE4'):
			        TOTAL_E.append(i.label)
			        n_elenode.append(4)
			        n_elements+=1
			        element_mark[i.label-1]=1
			#frame cycle
			for i_frame in range(0,len(odb.steps[stepname].frames)):
				
				#Detect whether the input frame is out of range
				try:
					TRY = odb.steps[stepname].frames[int(i_frame)]
				except:
					print "input frame exceeds the range of frames" 
					os._exit(0)
					break
				
				#Access a frame
				N_Frame = odb.steps[stepname].frames[int(i_frame)]
				print "Frame:",i_frame
				
				#create array for store result data temporarily
				# Vector-U,A,V,RF 
				L0=[] 
				# Tensors-S
				L1=[]
				# Tensors-LE
				L2=[]
				# Tensors-PE
				L3=[]
				# Scalars-PEEQ
				L4=[]
				for i in range(MLN): 
					L0.append([0,0,0])
					L1.append([0,0])
					L2.append([0,0])
					L3.append([0,0])
					L4.append([0,0])

				time1 = time()
				
				

				#Access Spatial displacement
				displacement = N_Frame.fieldOutputs['U']
				fieldValues = displacement.values
				for valueX in fieldValues :
					i = valueX.nodeLabel
					L0[i-1][0] = valueX.data[0]
					L0[i-1][1] = valueX.data[1]
					L0[i-1][2] = 0.0
					
				#SDV6
				SDV6 = N_Frame.fieldOutputs['SDV6']
				node_SDV6 = SDV6.getSubset(position=ELEMENT_NODAL)
				fieldValues = node_SDV6.values
				for valueX in fieldValues :
					L1[valueX.nodeLabel-1][0] += 1
					L1[valueX.nodeLabel-1][1] += valueX.data


				#SDV7
				SDV7 = N_Frame.fieldOutputs['SDV7']
				node_SDV7 = SDV7.getSubset(position=ELEMENT_NODAL)
				fieldValues = node_SDV7.values
				for valueX in fieldValues :
					L2[valueX.nodeLabel-1][0] += 1
					L2[valueX.nodeLabel-1][1] += valueX.data

				
				#SDV8
				SDV8 = N_Frame.fieldOutputs['SDV8']
				node_SDV8 = SDV8.getSubset(position=ELEMENT_NODAL)
				fieldValues = node_SDV8.values
				for valueX in fieldValues :
					L3[valueX.nodeLabel-1][0] += 1
					L3[valueX.nodeLabel-1][1] += valueX.data
				
					
				#SDV15
				SDV15 = N_Frame.fieldOutputs['SDV15']
				node_SDV15 = SDV15.getSubset(position=ELEMENT_NODAL)
				fieldValues = node_SDV15.values
				for valueX in fieldValues :
					L4[valueX.nodeLabel-1][0] += 1
					L4[valueX.nodeLabel-1][1] += valueX.data
	

				
				
				'''============================================================'''
				
				print "Partitionning model and writing vtk files ......"

				time1 = time()
				print "frame:",i_frame
				#Reorganization
				#Control&Storage
				#estimate whether the node has already existed
				stg_p = []
				#store the reorganized node for element
				stg_e = []
				#store the reorganized node for node
				stg_n = []
				for i in range(MLN):
					stg_p.append(-1)
				nodecount = 0
				#reorganize the node and element (reconstruct the mesh)

				M = range(0,n_elements)

				for i in M:
				    
					for j in range(n_elenode[i]):
						k = element[int(L_E[int(TOTAL_E[i]-1)])].connectivity[j] - 1
						if(stg_p[k] < 0): 
							stg_p[k] = nodecount
							stg_n.append(L[k]) 
							stg_e.append(nodecount)
							nodecount += 1
						else:
							stg_e.append(stg_p[k])
				#compute point quantity
				n_reop = len(stg_n)
				reop_N = range(0,len(stg_n))
			

				#create and open a VTK(.vtu) files

				outfile = open (os.path.join(vtk_path,odbname)+'_'+stepname+'_'+instancename+'f%03d'%int(i_frame)+'.vtu','w')
				
				#<VTKFile>, including the type of mesh, version, and byte_order
				outfile.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'+'\n')
				#<UnstructuredGrid>
				outfile.write('<UnstructuredGrid>'+'\n')
				#<Piece>, including the number of points and cells

				outfile.write('<Piece NumberOfPoints="'+str(n_reop)+'"'+' '+'NumberOfCells="'+str(n_elements)+'">'+'\n')

					
				print "Writing Nodes ......"
				#<Points> Write nodes into vtk files
				displacement = N_Frame.fieldOutputs['U']
				fieldValues = displacement.values
				outfile.write('<Points>'+'\n')
				outfile.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">'+'\n')
				for i in reop_N:
					nt = stg_n[i]
					k = node[stg_n[i]].label-1
					X,Y,Z = node[nt].coordinates[0]+L0[k][0]*0.5,node[nt].coordinates[1]+L0[k][1]*0.5,node[nt].coordinates[2]+L0[k][2]*0.5
					outfile.write(' '+'%11.8e'%X+'  '+'%11.8e'%Y+'  '+'%11.8e'%Z+'\n')			
				outfile.write('</DataArray>'+'\n')
				outfile.write('</Points>'+'\n')
				#</Points>

				print "Writing Results data ......"
				#<PointData> Write results data into vtk files
				outfile.write("<"+"PointData"+" "+"Scalars="+'"'+"Order Parameter,Stress_xx,Stres_yy,Stress_xy"+'"'+">"+'\n')
				#SDV15, <DataArray>
				outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Order Parameter"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
				for i in reop_N:
					k = node[stg_n[i]].label-1
					X = L4[k][1]/L4[k][0]
					outfile.write('%11.8e'%X+'\n')
				outfile.write('</DataArray>'+'\n')
				#SDV6, <DataArray>
				outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_xx"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
				for i in reop_N:
					k = node[stg_n[i]].label-1
					X = L1[k][1]/L1[k][0]
					outfile.write('%11.8e'%X+'\n')
				outfile.write('</DataArray>'+'\n')
				#SDV7, <DataArray>
				outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_yy"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
				for i in reop_N:
					k = node[stg_n[i]].label-1
					X = L2[k][1]/L2[k][0]
					outfile.write('%11.8e'%X+'\n')
				outfile.write('</DataArray>'+'\n')
				#SDV8, <DataArray>
				outfile.write("<"+"DataArray"+" "+"type="+'"'+"Float32"+'"'+" "+"Name="+'"'+"Stress_xy"+'"'+" "+"format="+'"'+"ascii"+'"'+">"+'\n')
				for i in reop_N:
					k = node[stg_n[i]].label-1
					X = L3[k][1]/L3[k][0]
					outfile.write('%11.8e'%X+'\n')
				outfile.write('</DataArray>'+'\n')
					#</DataArray>
			
				outfile.write("</PointData>"+'\n')
				#</PointData>
				
					
				print "Writing Cells ......"
				#<Cells> Write cells into vtk files
				outfile.write('<Cells>'+'\n')
				#Connectivity
				outfile.write('<DataArray type="Int32" Name="connectivity" format="ascii">'+'\n')
				offset=[0]
				for i in M:
				    if (n_elenode[i] == 3):
				        offset.append(offset[i]+3)
					for j in range(offset[i],offset[i+1]):
						outfile.write(str(stg_e[j])+' ')
				        outfile.write('\n')
				    elif (n_elenode[i] == 4):
				        offset.append(offset[i]+4)
					for j in range(offset[i],offset[i+1]):
						outfile.write(str(stg_e[j])+' ')
				        outfile.write('\n')	
				outfile.write('</DataArray>'+'\n')
				#Offsets
				outfile.write('<DataArray type="Int32" Name="offsets" format="ascii">'+'\n')
				for i in M:
				    outfile.write(str(offset[i+1])+'\n')
				outfile.write('</DataArray>'+'\n')
				#Type
				outfile.write('<DataArray type="UInt8" Name="types" format="ascii">'+'\n')
				for i in M:
				    if (n_elenode[i] == 3):
				        outfile.write('5\n')
				    elif (n_elenode[i] == 4):
				        outfile.write('9\n')
				outfile.write('</DataArray>'+'\n')
				outfile.write('</Cells>'+'\n')
				#</Cells>
				
				#</Piece>
				outfile.write('</Piece>'+'\n')
				#</UnstructuredGrid>
				outfile.write('</UnstructuredGrid>'+'\n')
				#</VTKFile>
				outfile.write('</VTKFile>'+'\n')
			
				outfile.close()

				
				'''====================================================================='''
				print "Creating .pvtu file for frame ", i_frame," ......"

			
			odb.close()

	print "Total time elapsed: ", time() - starttime, "s"
path=os.getcwd()
ConvertOdb2VtkP(path,'SingleCrack_UEL',path+'/Result/','0','100','1','0')	

		

