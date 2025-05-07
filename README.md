conda create -n OPLSCM5 python=3.11
conda activate OPLSCM5
have orca installed.
conda install -c conda-forge rdkit openbabel
pip install -e .
itp_rewrite -i carb -c CM5_charges.csv -o carb_orca.log
