# OPLSCM5
This package enables generation of GROMACS-compatible topology file for the OPLS-AA force field, using CM5 partial charges, denoted as OPLS-CM5 force field. 
input : pdb file
output : itp file

## Reference
This software utilizes a few archived or deprticated packages,particularly mktop and Orca2CM5charges, which needed significant modification. If you are using CNCstruc software, please ensure you cite the relevant papers for the aformentioned packages.
  1. LigParGen web server: an automatic OPLS-AA parameter generator for organic ligands  
     Leela S. Dodda  Israel Cabeza de Vaca  Julian Tirado-Rives William L. Jorgensen 
     Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W331–W336

  2. Charge Model 5: An Extension of Hirshfeld Population Analysis for the Accurate Description of Molecular Interactions in Gaseous and Condensed Phases
     Marenich, A. V.; Jerome, S. V.; Cramer, C. J.; Truhlar, D. G. J. Chem. Theory Comput. 2012, 8, 527– 541

  3. Evaluation of CM5 Charges for Condensed-Phase Modeling
     Jonah Z. Vilseck, Julian Tirado-Rives, and William L. Jorgensen J. Chem. Theory Comput., 2014, 10 (7), pp 2802–2812

  4. Evaluation of CM5 Charges for Nonaqueous Condensed-Phase Modeling
     Leela S. Dodda, Jonah Z. Vilseck, Kara J. Cutrona, and William L. Jorgensen J. Chem. Theory Comput., 2015, 11 (9), pp 4273–428
     
  5. MKTOP: a program for automatic construction of molecular topologies
     Ribeiro, A.A.S.T.; Horta, B.A.C.; de Alencastro, R.B.  J. Braz. Chem. Soc., 2008, 19 (7), pp 1433-1435, 

## Installation
Cloning the repository and creating a conda environment.

### Dependencies:
1 - ORCA software for partial charge calculation

2 - GROMACS, to detect the bonded interaction parameters from the GROMACS library.

### Installing the OPLSCM5 package
conda create -n OPLSCM5 python=3.11

conda activate OPLSCM5

conda install -c conda-forge rdkit openbabel

pip install -e .

## Application
Generating the topolgy file for the molecule.pdb:

OPLSCM5_gen -i molecule.pdb

