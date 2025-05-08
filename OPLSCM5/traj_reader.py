## this file includes multiple functions
import numpy as np
import os 
import  time
import subprocess
import linecache
import re
import pandas as pd
####### READERS
# this function reads and stores the coordinates, and outputs the DataFrame of coordinates and the box size
def pdb_reader (traj_file):
    atoms=[]
    f=open(traj_file,'r')
    output_order = ['residue_number', 'residue_name', 'atom_name', 'atom_number', 'x', 'y', 'z']
    for j,line in enumerate(f):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_info = {}
            # atom_info['record_name'] = line[0:6].strip()
            atom_info['atom_number'] = int(line[6:11].strip())
            atom_info['atom_name'] = line[12:16].strip()
            # atom_info['alternate_location'] = line[16].strip()
            atom_info['residue_name'] = line[17:20].strip()
            # atom_info['chain_id'] = line[21].strip()
            atom_info['residue_number'] = int(line[22:26].strip())
            # atom_info['icode'] = line[26].strip()
            atom_info['x'] = 0.1 * float(line[30:38].strip())
            atom_info['y'] = 0.1 * float(line[38:46].strip())
            atom_info['z'] = 0.1 * float(line[46:54].strip())
            # atom_info['occupancy'] = float(line[54:60].strip())
            # atom_info['temperature_factor'] = float(line[60:66].strip())
            # atom_info['element'] = line[76:78].strip()
            # atom_info['charge'] = line[78:80].strip()
            atoms.append({key : atom_info[key] for key in output_order})
    # data=num_conv_pd(np.array(coord),atom_type) ## save as dataframe
    return pd.DataFrame(atoms)
def gro_reader (traj_file,frame=1):
    f=open(traj_file,'r')
    lines = f.readlines()
    # coord=[]
    atoms = []
    # atom_type=[]
    atom_nums=atom_count(traj_file)
    frame_line=atom_nums+3       
    for line in lines[2:2+atom_nums]:
        atom_info = {}
        atom_info['residue_number'] = int(line[:5].strip())
        atom_info['residue_name'] = line[5:10].strip()
        atom_info['atom_name'] = line[10:15].strip()
        atom_info['atom_number'] = int(line[15:20].strip())
        atom_info['x'] = float(line[20:28].strip())
        atom_info['y'] = float(line[28:36].strip())
        atom_info['z'] = float(line[36:44].strip())
        if len(line) > 45:  # If velocities are present
            atom_info['vx'] = float(line[44:52].strip())
            atom_info['vy'] = float(line[52:60].strip())
            atom_info['vz'] = float(line[60:68].strip())
        atoms.append(atom_info)
    return pd.DataFrame(atoms)
def ndx_reader(data_pd,index_file,sel_ind):
    f1=open(index_file)
    atom_num=[]
    index_num=-1  ## number of index
    flag=0
    read_lines=''
    for j,line in enumerate(f1):
        if line.split()[0][0] == '[': 
            index_num=index_num+1
            if index_num==sel_ind: ## starts reading when we get to the desired index
                flag=1
                continue
            if index_num==sel_ind+1: ## stops reading when we pass the desired index
                break
        if flag:
            read_lines=read_lines+line
    atom_num=read_lines.split()
    atom_ind=[int(x) for x in atom_num]
    data_sel=data_pd.loc[atom_ind]
    return data_sel
## this function reads the xvg files, it neads the the filename and the columns as inputs. Outputs them as x and y. "xvg" reader
def ndx_writer(index_file,Data,group):
    f=open(index_file,"a")
    f.write('\n['+group+']\n')
    for iter,data_point in enumerate(Data):
        if (iter+1)%15==0:
            f.write('%5d \n' % data_point)
        else:
            f.write('%5d      ' % data_point)
    f.close()
def xvg_reader(filename,indice):
      spec_char=['@','#','&']
      x=[]
      y=[]
      with open(filename) as f:
            for j,line in enumerate(f):
                  if line.split()[0][0]  in spec_char:
                        continue
                  x.append(float(line.split()[indice[0]-1]))
                  y.append(float(line.split()[indice[1]-1]))
      return x,y
####### FUNCTIONS
## this function counts total number of particles in the box
def atom_count(traj_file):
    if traj_file[-3:]=='gro':
        f1=open(traj_file)
        for j,line in enumerate(f1):
            if j==1:
                atom_nums=int(line)
                break
        return atom_nums
def frame_count(traj_file):
    atom_nums=atom_count(traj_file)
    line_num = sum(1 for line in open(traj_file))  ## for very big files , we recommed linux users to uncomment the above
    frame_lines=atom_nums+3
    return int(line_num/frame_lines)    

## this function converts coordinate 
def num_conv_pd(coordinate,atom_type):
    dicts={'Bead': atom_type,'X': coordinate[:,0], 'Y':  coordinate[:,1], 'Z':  coordinate[:,2]}
    data=pd.DataFrame(dicts)
    data.index+=1
    return data
## this function writes in the gro format
def gro_writer_old(traj_file,CG_data,box_size,res_name,frame=1):
    f=open(traj_file,"a")
    f.write('frame number is %d\n' % frame)
    f.write(' %6d\n' % len(CG_data))

    for i in range(len(CG_data)):
        f.write('%5i%-5s%5s%5i%8.3f%8.3f%8.3f\n' % (1,res_name,CG_data.loc[i+1,'Bead'],float(i+1),CG_data.loc[i+1,'X'],CG_data.loc[i+1,'Y'],CG_data.loc[i+1,'Z']))
    f.write('%10s%20s%20s\n' % (box_size.split()[0],box_size.split()[1],box_size.split()[2])) 
    f.close()
    return f
def pdb_writer(traj_file,data):
    f=open(traj_file,"a")
    for i in data.index.tolist():
        atom_line = "{:<6}{:>5} {:<4}{:<1}{:<3} {:<1}{:>4}{:<1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2}{:>6.2}          {:>2}{:>2}".format(
        data.loc[i,'record_name'] if 'record_name' in data.columns else 'ATOM',
        data.loc[i,'atom_number'],
        data.loc[i,'atom_name'],
        data.loc[i,'alternate_location'] if 'alternate_location' in data.columns else ' ',
        data.loc[i,'residue_name'],
        data.loc[i,'chain_id']  if 'chain_id' in data.columns else ' ',
        data.loc[i,'chain_number'] ,
        data.loc[i,'icode'] if 'icode' in data.columns else ' ',
        data.loc[i,'x'],
        data.loc[i,'y'],
        data.loc[i,'z'],
        data.loc[i,'occupancy'] if 'occupancy' in data.columns else ' ' ,
        data.loc[i,'temperature_factor'] if 'temperature_factor' in data.columns else ' ' ,
        data.loc[i,'element'] if 'element' in data.columns else ' ' ,
        data.loc[i,'charge'] if 'charge' in data.columns else ' ')
        f.write(atom_line + "\n")

    f.close()
    return f

def gro_writer(traj_file,data,box_dimensions=[10,10,10]):
    # res_name="WAT"
    # atom_name="OW"
    # f=open(traj_file,"a")
    f=open(traj_file,"w")
    f.write("Generated by script\n")
    f.write("{:5d}\n".format(len(data)))
    for i in data.index.tolist():
        atom_line = "{:5d}{:<5}{:>5}{:5d}{:8.3f}{:8.3f}{:8.3f}".format(
        data.loc[i,'residue_number'] ,
        data.loc[i,'residue_name'],
        data.loc[i,'atom_name'],
        data.loc[i,'atom_number'],
        data.loc[i,'x'],
        data.loc[i,'y'],
        data.loc[i,'z'])
        f.write(atom_line + "\n")
    f.write("{:10.5f}{:10.5f}{:10.5f}\n".format(*box_dimensions))
    return f