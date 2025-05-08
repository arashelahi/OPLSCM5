import numpy as np
import pandas as pd
def itp_reader(file): ## draw the opls types of atoms
    flag=0
    opls_types=[]
    f=open(file,'r')

    for line in f:
        if line.split()[:3]==['[','atoms',']']:
            flag=1
            continue
        if line.split()[:3]==['[','bonds',']']:
            break
        if ((line.split()==[]) or (line.split()[0]==';')):
            continue
        if flag==1:
            opls_types.append(line.split()[1])
    return opls_types

def At_type_finder(nonbond_file,itp_file): ## make a DataFrame of atom numbres, opls names and atomtypes
    opls_types=itp_reader(itp_file)
    opls_types=opls_types
    atom_types=[]
    for at_ind in range(len(opls_types)):
        f=open(nonbond_file,'r')
        for line in f:
            if line.split()[0]==opls_types[at_ind]:
                atom_types.append(line.split()[1])
                continue
        f.close()
    Dicts={'At_numbers' : [i+1 for i in range(len(opls_types))],'opls_types' : opls_types,'At_types' : atom_types}
    datas=pd.DataFrame(Dicts)
    datas.index=datas.index+1
    return datas

def atom_itp(nonbond_file,itp_old_file,itp_new_file,Charges,mol_name,nrexcl=3): ## creates the Mol and atom section in the itp file
    itp_old_read=open(itp_old_file,'r') ## reading the old itp file
    itp_new_read=open(itp_new_file,'a') ## reading writtable new itp file
    datas=At_type_finder(nonbond_file,itp_old_file)
    flag=0
    iter=0

    for line in itp_old_read:
        text=line
        if line.split()[:3]==['[','moleculetype',']']:  
            flag='molectype' 
            text=line+";%-6s%10s\n" % ('Name','nrexcl')
        if line.split()[:3]==['[','atoms',']']:
            flag='atomtype'
            text=line+";  %6s%10s%10s%10s%10s%15s%20s%15s\n" % ('nr','type','resnr','residue','atom','cgnr','charge','mass')
        if line.split()[:3]==['[','bonds',']']:
            itp_new_read.close()
            break
        if flag==0:  continue
        if flag=='molectype':
            if ((line.split()==[]) or (line.split()[0]==';')):
                continue
            elif line.split()[:3]!=['[','moleculetype',']']:
                text='%-8s%2s\n' % (mol_name,nrexcl)
        if flag=='atomtype':
            if ((line.split()==[]) or (line.split()[0]==';')):
                continue
            elif line.split()[:3]!=['[','atoms',']']:
                iter=iter+1
                text="   %6s%10s%10s%10s%10s%15s%20s%15s\n" % (datas.loc[iter,'At_numbers'],datas.loc[iter,'opls_types'],\
                    line.split()[2],mol_name,datas.loc[iter,'At_types'],datas.loc[iter,'At_numbers'],Charges[iter-1],line.split()[7])
        itp_new_read.write(text)

def bond_itp(nonbond_file,bond_file,itp_old_file,itp_new_file): ## creates the bond section for the itp file
    itp_old_read=open(itp_old_file,'r')
    itp_new_read=open(itp_new_file,'a')
    datas=At_type_finder(nonbond_file,itp_old_file)
    flag=0
    
    for line in itp_old_read:
        text=line
        if line.split()[:3]==['[','bonds',']']: 
            flag='bondtype'
            text=line+"; %4s%6s%6s%15s%15s\n" % ('i','j','func','b','kb')
        if flag==0:continue
        elif line.split()[:3]==['[','angles',']']:
            itp_new_read.close()
            break
        if flag=='bondtype':
            if ((line.split()==[]) or (line.split()[0]==';')):
                continue
            elif line.split()[:3]!=['[','bonds',']']:
                atom1=datas.loc[int(line.split()[0]),'At_types']
                atom2=datas.loc[int(line.split()[1]),'At_types']
                flag_b=0
                f=open(bond_file,'r')
                for line_b in f:
                    if line_b.split()[:3]==['[','bondtypes',']']:
                        flag_b='bondtype'
                        continue
                    if flag_b=='bondtype':
                        if ((line_b.split()==[]) or (line_b.split()[0]==';')):continue
                        if line_b.split()[0]=='[':break
                        if (((line_b.split()[0]==atom1) and (line_b.split()[1]==atom2)) or\
                            ((line_b.split()[0]==atom2) and (line_b.split()[1]==atom1))):
                            text="%4s%6s%6s%15s%15s\n" % (line.split()[0],line.split()[1]\
                                                    ,line_b.split()[2],line_b.split()[3],line_b.split()[4] )
                            f.close()
                            break           
        itp_new_read.write(text)

def angle_itp(nonbond_file,bond_file,itp_old_file,itp_new_file):
    itp_old_read=open(itp_old_file,'r')
    itp_new_read=open(itp_new_file,'a')
    datas=At_type_finder(nonbond_file,itp_old_file)
    flag=0
    for line in itp_old_read:
        text='line'
        if line.split()[:3]==['[','angles',']']:
            text=line+"; %4s%6s%6s%6s%15s%15s\n" % ('i','j','k','func','b','kb')         
            flag='angletype'
        if line.split()[:3]==['[','dihedrals',']']:
                itp_new_read.close()
                break
        if flag==0:
            continue
        if flag=='angletype':
            if ((line.split()==[]) or (line.split()[0]==';')):continue
            elif line.split()[:3]!=['[','angles',']']:
                atom1=datas.loc[int(line.split()[0]),'At_types']
                atom2=datas.loc[int(line.split()[1]),'At_types']
                atom3=datas.loc[int(line.split()[2]),'At_types']
                flag_b=0
                f=open(bond_file,'r')
                find_stat=0
                for line_b in f:
                    if line_b.split()[:3]==['[','angletypes',']']:
                        flag_b='angletype'
                        continue
                    if flag_b=='angletype':
                        if ((line_b.split()==[]) or (line_b.split()[0]==';')):
                            continue
                        if line_b.split()[0]=='[':break
                        if (((line_b.split()[0]==atom1) and (line_b.split()[1]==atom2) and (line_b.split()[2]==atom3)) or\
                            ((line_b.split()[0]==atom3) and (line_b.split()[1]==atom2) and (line_b.split()[2]==atom1))):
                            find_stat=1
                            text="%4s%6s%6s%6s%15s%15s\n" % (line.split()[0],line.split()[1],line.split()[2]\
                                                    ,line_b.split()[3],line_b.split()[4],line_b.split()[5] )
                            f.close()
                            break
                if find_stat==0:
                    text="%4s%6s%6s%6s%15s%15s%10s\n" % (line.split()[0],line.split()[1],line.split()[2],\
                        atom1,atom2,atom3,'NOT FOUND' )
                    f.close()
                    break
        itp_new_read.write(text)

def dih_itp(nonbond_file,bond_file,itp_old_file,itp_new_file,dihedral_type='proper'):
    itp_old_read=open(itp_old_file,'r')
    itp_new_read=open(itp_new_file,'a')
    datas=At_type_finder(nonbond_file,itp_old_file)
    flag=0
    dih_count=0

    for line in itp_old_read:
            
        text=line
        if line.split()[:3]==['[','dihedrals',']']:
            text=line+"%4s%6s%6s%6s%6s%15s%15s%15s%15s%15s%15s\n" % (';ai','aj','ak','al','funct','c0','c1','c2','c3','c4','c5')
            flag='dihedtype'
        if flag==0:continue
        if line.split()[:3]==['[','pairs',']']: ## the section after dihedral is not always same
            itp_new_read.close()
            break
        if line.count('improper_')==1:
            flag='dihedtype_imp'
        flag_b=0
        if flag=='dihedtype_imp':
            if ((line.split()==[]) or (line.split()[0]==';')):continue
            f=open(bond_file,'r')
            for line_b in f:
                if line_b.split()[:3]==['[','dihedraltypes',']']:
                    dih_count=dih_count+1
                    if dih_count==1:
                        continue
                    flag_b='dihedtype_imp'
                    continue
                if flag_b=='dihedtype_imp':
                    if ((line_b.split()==[]) or (line_b.split()[0]==';')):
                        continue
                    if line_b.split()[0]=='[':
                        break
                    if line.split()[5]==line_b.split()[1]:
                        text="%4s%6s%6s%6s%6s%15s%15s%15s\n" % (line.split()[0],line.split()[1],line.split()[2],line.split()[3]\
                            ,str(4),line_b.split()[2],line_b.split()[3],line_b.split()[4])
                        f.close()
                        break
        # if flag=='dihedtype_imp':
        #     itp_new_read.write(text)
        #     continue
        if flag=='dihedtype':
            if ((line.split()==[]) or (line.split()[0]==';')):continue
            elif line.split()[:3]!=['[','dihedrals',']']:
                atom1=datas.loc[int(line.split()[0]),'At_types']
                atom2=datas.loc[int(line.split()[1]),'At_types']
                atom3=datas.loc[int(line.split()[2]),'At_types']
                atom4=datas.loc[int(line.split()[3]),'At_types']
                Torsion_atoms=[atom1,atom2,atom3,atom4]
                flag_b=0
                f=open(bond_file,'r')
                find_stat=0
                for line_b in f:
                    if line_b.split()[:3]==['[','dihedraltypes',']']:
                        flag_b='dihedtype'
                        continue
                    if flag_b=='dihedtype':
                        if ((line_b.split()==[]) or (line_b.split()[0]==';')):
                            continue
                        if line_b.split()[0]=='[':
                            break
                        if (((line_b.split()[0]==atom1) and (line_b.split()[1]==atom2) and (line_b.split()[2]==atom3) and (line_b.split()[3]==atom4)) or\
                            ((line_b.split()[0]==atom4) and (line_b.split()[1]==atom3) and (line_b.split()[2]==atom2) and (line_b.split()[3]==atom1))):
                            find_stat=1
                            text="%4s%6s%6s%6s%6s%15s%15s%15s%15s%15s%15s\n" % (line.split()[0],line.split()[1],line.split()[2],line.split()[3]\
                                ,line_b.split()[4],line_b.split()[5],line_b.split()[6],line_b.split()[7],line_b.split()[8],line_b.split()[9],line_b.split()[10])
                            f.close()
                            break

                if find_stat==0:
                    flag_b=0
                    f=open(bond_file,'r')
                    find_stat=0
                    for line_b in f:
                        if line_b.split()[:3]==['[','dihedraltypes',']']:
                            flag_b='dihedtype'
                            continue
                        if flag_b=='dihedtype':
                            if ((line_b.split()==[]) or (line_b.split()[0]==';')):
                                continue
                            if line_b.split()[0]=='[':
                                break
                            if (np.sum(np.array(line_b.split()[:4])==np.array(Torsion_atoms))==3 or\
                                np.sum(np.array(line_b.split()[:4])==np.array(Torsion_atoms[::-1]))==3):
                                if (np.sum(np.array([x[0] for x in line_b.split()[:4]])==np.array([x[0] for x in Torsion_atoms]))==4 or\
                                    np.sum(np.array([x[0] for x in line_b.split()[:4]])==np.array([x[0] for x in Torsion_atoms[::-1]]))==4):

                                    if Torsion_atoms[np.argmax(np.array(line_b.split()[:4])!=np.array(Torsion_atoms))] in ['CO','C']:
                                        Torsion_atoms[np.argmax(np.array(line_b.split()[:4])!=np.array(Torsion_atoms))]='CT'
                                    elif Torsion_atoms[::-1][np.argmax(np.array(line_b.split()[:4])!=np.array(Torsion_atoms[::-1]))] in ['CO','C']:
                                        Torsion_atoms[np.argmax(np.array(line_b.split()[:4][::-1])!=np.array(Torsion_atoms))]='CT'

                            if (np.sum(np.array(line_b.split()[:4])==np.array(Torsion_atoms))==4 or\
                                np.sum(np.array(line_b.split()[:4])==np.array(Torsion_atoms[::-1]))==4):
                                find_stat=1
                                text="%4s%6s%6s%6s%6s%15s%15s%15s%15s%15s%15s\n" % (line.split()[0],line.split()[1],line.split()[2],line.split()[3]\
                                ,line_b.split()[4],line_b.split()[5],line_b.split()[6],line_b.split()[7],line_b.split()[8],line_b.split()[9],line_b.split()[10])
                                f.close()
                                break

                if find_stat==0:
                    line_b='pass pass pass pass 3        0.00000        0.00000        0.00000        0.00000        0.00000        0.00000'
                    text="%4s%6s%6s%6s%6s%15s%15s%15s%15s%15s%15s\n" % (line.split()[0],line.split()[1],line.split()[2],line.split()[3]\
                    ,line_b.split()[4],line_b.split()[5],line_b.split()[6],line_b.split()[7],line_b.split()[8],line_b.split()[9],line_b.split()[10])
                    # text="%4s%6s%6s%6s%6s%15s%15s%15s%15s\n" % (line.split()[0],line.split()[1],line.split()[2],line.split()[3],\
                    #     atom1,atom2,atom3,atom4,'NOT FOUND' )
                    f.close()
        itp_new_read.write(text)

def rest_itp(itp_old_file,itp_new_file):
    itp_old_read=open(itp_old_file,'r')
    itp_new_read=open(itp_new_file,'a')
    flag=0
    for line in itp_old_read:
        text=line
        if line.split()[:3]==['[','pairs',']']:
            # itp_new_read.(line)
            flag='pairtype'     
        if flag==0: continue   
        if line.split()[:3]==['[','system',']']:
            itp_new_read.close()
            break
        itp_new_read.write(text)


