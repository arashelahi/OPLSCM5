from CNCstruc.structure import CNC_class as CNC
from CNCstruc.utils import traj_reader as trj
# from CNCstruc.analysis import Indexing as indx_gen
import pandas as pd
import os 
import numpy as np
from openbabel import pybel
# from rdkit import Chem
# from rdkit.Chem import AllChem


def molecule_conversion(parent_molecule , reference_atom , atoms_to_delete , resid_range , func_group_init):
    '''
    This function converts the reactive groups to the modifying functional group:
    input: 
    parent_molecule : the datafile for the parent molecule
    reference_atom : the reference atom where the functional group and the reactive groups meet
    atoms_to_delete : the atom names which should be deleted from the functional group.
    resid_range : the residues where the functional groups will be applied upon
    func_group_init : the functional group dataframe
    # Outputs:
    the resulting data file
    '''
    for resid in resid_range:
        # The translational movement of the functional group replacement
        relative_position = parent_molecule[(parent_molecule['residue_number'] == resid) & 
                                            (parent_molecule['atom_name'] == reference_atom)].loc[:,'x':'z'].values - func_group_init[func_group_init['atom_name'] == reference_atom].loc[:,'x':'z'].values
        # remove the merging atoms from the functional group
        func_group = func_group_init[func_group_init['atom_name'] != reference_atom]
        # the movement of the functional group
        func_group.loc[:,'x':'z'] = func_group_init.loc[:,'x':'z'] + relative_position
        # making the functional group properties aligned with the parent molecule
        func_group[['residue_number', 'residue_name' , 'chain_number']] = parent_molecule[parent_molecule['residue_number'] == resid][['residue_number', 'residue_name','chain_number']].reset_index(drop=True)
        # delete the atoms which are converted
        data_to_delete = parent_molecule[(parent_molecule['residue_number'] == resid) & parent_molecule['atom_name'].isin(atoms_to_delete)].index
        parent_molecule = parent_molecule.drop(data_to_delete) # delete the atoms of the alcohol groups
        # The reordering is based on atom number. Shift the atom numbers of partciles after the removal as well as the functional group
        parent_molecule.loc[data_to_delete[0] + 1 :,'atom_number'] = np.arange(1 , len(parent_molecule.loc[data_to_delete[0] + 1 :,'atom_number'])+1 )\
                                                                        + parent_molecule.loc[data_to_delete[0] - 1,'atom_number'] + len(func_group)  # then the next atom will be shifted by the number of replacing atoms
        func_group.loc[:,'atom_number'] = np.arange(1 , len(func_group)+1 ) + parent_molecule.loc[data_to_delete[0] - 1 ,'atom_number']
        # Concat the functional groups to the cellulose
        parent_molecule = pd.concat([parent_molecule , func_group])
        # Reorder based on the atom number
        parent_molecule = parent_molecule.sort_values(by = 'atom_number' , ascending = True).reset_index(drop = True)
    return parent_molecule

def functionalization (parent_molecule , func , side):
    # The functional group data and the residues to replace
    func_group_init , resid_range = _data_preparation(parent_molecule , func , side)
    # The atom to be replaced and the atoms to be deleted
    reference_atom = 'C6'
    if func == 'COOH':
        reference_atom = 'C6'
        atoms_to_delete = ('H61', 'H62', 'O6' , 'HO6')
    else:
        # if the group is alkyl, it has to be applied on the carboxylated cellulose
        parent_molecule = functionalization (parent_molecule , 'COOH' , side)
        atoms_to_delete = ('O','H')
    functionalized_CNC = molecule_conversion(parent_molecule , reference_atom , atoms_to_delete , resid_range , func_group_init)
    return functionalized_CNC

def _data_preparation(parent_molecule , func , side):
    functional_file = f'./data/input/Groups/{func}_{side}.pdb'
    functional_data = trj.pdb_reader(functional_file)
    if side == 'right':
        resid_range = np.arange(parent_molecule['residue_number'].min() , parent_molecule['residue_number'].max() , 2)
    else:
        resid_range = np.arange(parent_molecule['residue_number'].min() + 1, parent_molecule['residue_number'].max()+1 , 2)
    return functional_data , resid_range
def CNC_creation (parent_molecule , func , left_molecules , right_molecules):
    # left_molecules = [first_elem[0] for first_elem in list(CNC_group.layers.values())[:-1]]
    # right_molecules = [first_elem[-1] for first_elem in list(CNC_group.layers.values())[0:]]
    new_CNC = pd.DataFrame(columns=parent_molecule.columns)
    for chain_num in np.unique(parent_molecule['chain_number']):
        current_data = parent_molecule[parent_molecule['chain_number'] == chain_num]
        if chain_num in left_molecules:
            side = 'left'
        elif chain_num in  right_molecules:
            side = 'right'
        else:
            new_CNC = pd.concat([new_CNC ,current_data ])
            continue
        current_data = functionalization(current_data , func, side)
        new_CNC = pd.concat([new_CNC ,current_data ])
    new_CNC['atom_number'] = new_CNC.reset_index(drop = True).index + 1
    return new_CNC.reset_index(drop = True)
def material_prep(CNC_group , func):
    ''' This function gets the CNC object and the functional group and creates the gro files for the surface-modified ones.
    inputs:
    CNC_group : The CNC class
    func : the functional groups which could be COOH, ethy, or butyl.
    '''
    Data = CNC_group.data
    cellulose_data = Data[Data['chain_number'] == 1]
    side_vec = ['left' , 'right']
    filepath =  f'./data/input/Modified/{func}/coordinate/'
    os.makedirs(filepath, exist_ok=True)
    for side in side_vec:
        cellulose_func = functionalization (cellulose_data , func , side)
        filename = filepath + f'cellulose_{func}_{side}'
        
        trj.gro_writer(f'{filename}.gro', cellulose_func)
        mol = next(pybel.readfile("gro",f'{filename}.gro'))
        mol.write("xyz", f"{filename}.xyz", overwrite=True)
        
    
    ### Making the CNC.
    left_molecules = [first_elem[0] for first_elem in list(CNC_group.layers.values())[:-1]]
    right_molecules = [first_elem[-1] for first_elem in list(CNC_group.layers.values())[0:]]
    CNC_COOH_Data =  CNC_creation (Data , func ,left_molecules , right_molecules)
    trj.gro_writer(f'{filepath}/CNC_{func}.gro' , CNC_COOH_Data)


