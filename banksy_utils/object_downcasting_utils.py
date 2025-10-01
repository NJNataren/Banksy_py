# Author: Nathalie Nataren
# Date: 30/09/
#  These small helper functions can be used to check the float size of the .obsm slot of AnnData objects and then downcast this data to reduce RAM usage

import numpy as np

# Function to check the float size of the adata.obsm data
def check_float(adata):
    for key in adata.obsm.keys():
        arr = adata.obsm[key]
    if isinstance(arr, np.ndarray):
        print(f"{key}: shape={arr.shape}, dtype={arr.dtype}")

# Function to downcast adata.obsm to 32 bit and check that the resulting array is an N-dimensional array
def downcast_float(adata, float_type):
    float_type = np.dtype(float_type)
    for key in adata.obsm.keys():
        arr = adata.obsm[key]
    if isinstance(arr, np.ndarray) and np.issubdtype(arr.dtype, np.floating): #check if the array is an N-dimensional array
        adata.obsm[key] = arr.astype(float_type)
        
        print(f"{adata} object downcast to {float_type}")  