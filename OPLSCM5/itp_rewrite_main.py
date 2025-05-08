import pandas as pd
import numpy as np
import os
import argparse
import subprocess
import OPLSCM5
import OPLSCM5.opls_rewrite as opls_itp_rewrite
import OPLSCM5.Orca2CM5charges as CM5

def get_gromacs_top_dir():
    
        # Run the shell command and capture output
    for cmd in [["gmx_mpi"] , ["gmx"]]:
    # cmd = ["gmx_mpi"]
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            for line in result.stderr.splitlines():
                if "Data prefix" in line:
                    # Split by colon and strip whitespace
                    prefix = line.split(":", 1)[1].strip()
                    return os.path.join(prefix, "share", "gromacs", "top" , "oplsaa.ff")

        except FileNotFoundError:
            print("Error: gmx_mpi not found.")
        except subprocess.CalledProcessError:
            print("Error: Failed to run gmx_mpi -help.")

def orca_inp_prep(molecfile):
    # Create the input file for ORCA

    text=f"!RKS RIJCOSX M062X cc-pVDZ PMODEL PAL6\n%maxcore 10000\n% output\n\n        Print[P_hirshfeld] 1\n\n end\n\n*xyzfile 0 1 {molecfile}.xyz\n\n"
    with open(f'{molecfile}_orca.inp', 'w') as f:
        f.write(text)
def main():
    # This is the main function that will be called when the script is run
# if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='usage of MKTOP and and ORCA')

    parser.add_argument("-i", "--molecfile", help="the top file from the mktop", type=str)
    args = parser.parse_args()
    mol_name=args.molecfile[:-4]
    itp_old_read='%s.top' % mol_name
    itp_new_read='%s.itp' % mol_name
    # charge_file=args.chargefile
    charge_file = '%s_charges.csv' % mol_name
    # ffpath = args.ffpath
    
    # Doing the Orca analysis
    print('####### Doing Orca calculation ################')
    print('creating the xyz file for orca')
    os.system(f'obabel {mol_name}.pdb -O {mol_name}.xyz')
    print('Running orca')
    orca_inp_prep(mol_name)
    ORCA_EXE = subprocess.check_output("which orca", shell=True, text=True).strip()
    os.system(f'{ORCA_EXE} {mol_name}_orca.inp > {mol_name}_orca.log')
    a0,rd,pt = CM5.LoadModel()
    data = CM5.GetLogFile(f'{mol_name}_orca.log',pt,rd)
    qcm5 = CM5.HirshfeldToCM5(data,a0)
    qcm5.to_csv(charge_file,index=False,float_format='%7.5f')
    # print(qcm5)

    print('####### Performing MKTOP analysis ################')
    GMX_DIR = get_gromacs_top_dir()
    # print(GMX_DIR)
    nonbond_file=os.path.join(GMX_DIR,'ffnonbonded.itp')
    bond_file=os.path.join(GMX_DIR,'ffbonded.itp')
    # print(bond_file)

    ## Creating the pdb file and top files using mktop and obabel
    # os.system(f'obabel {mol_name}.xyz -O {mol_name}.pdb')
    mktop_path = os.path.join(os.path.dirname(OPLSCM5.__file__), "mktop.pl")
    os.system(f'{mktop_path} -i {mol_name}.pdb -o {mol_name}.top -ff opls -conect no')

    if os.path.exists(str(itp_new_read)):
        os.remove(itp_new_read) 
    print('########## creating the new itp file ################')
    CM5_data=pd.read_csv(charge_file)
    Charges=np.array(CM5_data['1.20*CM5'])
    Charges[-1]=Charges[-1]-np.sum(Charges)
    Charges[-1]=float("%7.5f" % Charges[-1])
    #Charges=np.zeros(195)
    Datas=opls_itp_rewrite.At_type_finder(nonbond_file,itp_old_read)
    opls_itp_rewrite.atom_itp(nonbond_file,itp_old_read,itp_new_read,Charges,mol_name)
    opls_itp_rewrite.bond_itp(nonbond_file,bond_file,itp_old_read,itp_new_read)
    opls_itp_rewrite.angle_itp(nonbond_file,bond_file,itp_old_read,itp_new_read)
    opls_itp_rewrite.dih_itp(nonbond_file,bond_file,itp_old_read,itp_new_read)
    # opls_itp_rewrite.dih_itp(nonbond_file,bond_file,itp_old_read,itp_new_read,'improper')
    opls_itp_rewrite.rest_itp(itp_old_read,itp_new_read)

if __name__ == "__main__":
    main()