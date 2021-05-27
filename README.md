# This repo contains all code developed to pull relevant structural features down from EmbA and EmbB, as well as all notebooks for optimising and testing the 4 machine learning models    
### Conceptually, the code for rifampicin (not uplaoded) and ethambutol is identical, however, ethambutol repo was made accessible due to tidiness and better annotations   
### No clinicla data in the project was made public, and therefore, the machine learning models do not run - however, the Random Forest notebook does contain cached outputs that were saved from
### actual runs, and thus model performances can be analysed     


## Notes:   
#### All files of interest are under the 'notebooks' directory   
#### Py files are identical to their ipynb counterparts, and can in fact be run in absence of clinical data if the correct software is downloaded
#### Structural_attributes pulls down distances, depth, and secondary strcuture from the pdb files
#### thermo_features pulls relevant data from foldx output files
#### Amino_acid_attributes calculates changes in biochemical metrics on mutation, and determines evolutionary conservation via MAPP
#### Build_ML_df.ipynb essentially merges the clinical and structural data, and outputs a dataframe in the correct format for training 
