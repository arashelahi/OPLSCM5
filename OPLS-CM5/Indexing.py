'''
This module is used for making the index files for any arbitrary choice of atoms. 

'''

# import traj_reader as trj
# import CNC_class as cnc
import os
import numpy as np
import matplotlib.pyplot as plt
import pathlib
from CNCstruc.utils import traj_reader as trj

def remove_file(file_path):
    if os.path.isfile(file_path):
        os.remove(file_path)
    
def ndx_making(CNC_group, feature , ndx_file = None , output_path = ''):
    # remove_file(ndx_file)
    if feature in ['glycosidic', 'alcohols', 'twist']:
        dihed_make_ndx(CNC_group , feature , output_path)
    elif feature == 'H_bonds':
        hb_make_ndx(CNC_group , feature , output_path)
    elif feature == 'unit_cell':
        for feature_type in CNC_group.descriptor[feature].keys():
            unit_cell_make_ndx(CNC_group , feature , feature_type , output_path)
        # self.unit_angle_make_ndx(feature)
    else:
        raise ValueError('Invalid feature name. Please choose from the following: glycosidic, alcohols, twist, H_bonds, or unit_cell.')

def dihed_make_ndx(CNC_group , feature , output_path):
    ''' 
    This function creates the index file for the dihedral angles. 
    The dihedral angles are extracted from the CNC_class object. 
    '''
    ndx_file = output_path + feature+'_index.ndx'
    remove_file(ndx_file)
    for feature_type in CNC_group.descriptor[feature].keys():
        if feature in ['twist' , 'glycosidic']:
            trj.ndx_writer(ndx_file,CNC_group.descriptor[feature][feature_type],feature_type)
        elif feature == 'alcohols':
            data_right = [x for iter_x , x in enumerate(CNC_group.descriptor['alcohols'][feature_type]) if iter_x % 8 >= 4 and iter_x % 8 < 8]
            data_left = [x for iter_x , x in enumerate(CNC_group.descriptor['alcohols'][feature_type]) if iter_x % 8 < 4]
            if CNC_group.ff=='Charmm':
                reserved_data = data_left.copy()
                data_left = data_right.copy()
                data_right = reserved_data.copy()
            trj.ndx_writer(ndx_file,data_right,f'{feature_type}_left')
            trj.ndx_writer(ndx_file,data_left,f'{feature_type}_right')


def unit_cell_make_ndx(CNC_group , feature , feature_type , output_path):
    ndx_file = output_path + f'unit_cell_{feature_type}.ndx'
    remove_file(ndx_file)
    for cell_elem in CNC_group.descriptor[feature][feature_type].keys():
        trj.ndx_writer(ndx_file,CNC_group.descriptor[feature][feature_type][cell_elem], cell_elem)


def hb_make_ndx(CNC_group , feature , output_path):
    for feature_type in CNC_group.descriptor[feature].keys():
        ndx_file = output_path + f'{feature_type}_index.ndx'
        remove_file(ndx_file)
        block_size = len(CNC_group.ATOM_TYPES[feature][feature_type])*len(CNC_group.resid_vec)
        if feature_type in ['O2H_O6','O3H_O5']:
            layer_vec = CNC_group.layer_vec
        else:
            layer_vec = CNC_group.layer_vec[1:-1]
        iter = 0
        for layer in layer_vec:
            if feature_type in ['O2H_O6','O3H_O5']:
                chain_number_vec = CNC_group.layers[layer][1:-1]
            else:
                chain_number_vec = CNC_group.layers[layer][1:-2]
            for chain_number_iter,chain_number in enumerate(chain_number_vec):
                data = CNC_group.descriptor[feature][feature_type][iter*block_size:(iter+1)*block_size]
                trj.ndx_writer(ndx_file,data,"%s_ch%d_%s" % (layer,chain_number_iter,feature_type)) if data else None
                iter += 1


