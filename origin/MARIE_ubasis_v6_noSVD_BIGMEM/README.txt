this verion of the code was modified to minimize the GPU memory used at any one time. This comes at some cost of 
computation time, especially in the solve (MARIE) section of the code.

Therefore, only use this version of the code if memory is an issue, e.g. if MARIE_ubasis_v4_noSVD crashes.