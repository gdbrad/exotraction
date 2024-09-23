import numpy as np 
import gvar as gv 
import matplotlib 
import h5py as h5 
import pathlib 
from typing import Literal,NamedTuple,List,Tuple,Iterable


def C_2pt_master(
        C_2pt: pathlib.Path ,
        n_skip: List[Tuple[int]],
        nvecs:List[Tuple[int]]):
    ''''
    Returns: 
    h5 file with dimensions
    N_t x t_srcs x N_cfg x Num(nvecs)
    96 x 96 x 200 x 4
    '''
    

def C_2pt_raw(
        master_corr: np.ndarray,
        n_skip: List[Tuple[int]],
        nvecs:List[Tuple[int]]):
        for i,nvec in enumerate(nvecs):
              
                raw_corr =  master_corr[i,:,:,:].squeeze()
                for n in n_skip:
                        corr = raw_corr[:,::n,:]
                        corr = corr.mean(axis=1) 
        
                




