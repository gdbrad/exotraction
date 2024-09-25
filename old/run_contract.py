import ingest_data
import sys
import os
import numpy as np 
from pathlib import Path
import contract
# from pion_computation import compute_pion
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import functools
import h5_utils

from gamma import gamma
from contract import contract_pion_t, contract_vector_t
# from ingest_data import load_elemental, load_peram, reverse_perambulator_time
from ingest_data import load_elemental, load_peram, reverse_perambulator_time

h5_path = os.path.abspath(os.path.join('..','..','exolaunch'))

# Add the directory to sys.path
if h5_path not in sys.path:
    sys.path.append(h5_path)
h5_path
# get these from ens yml file 
Lt=96
nvec=32

for peram_file in os.listdir(os.path.join(h5_path,'perams_sdb',f'numvec{nvec}')):
    cfg_id = peram_file.strip('perams_')
    print(cfg_id)
    for meson_file in os.listdir(os.path.join(h5_path,'meson_sdb',f'numvec{nvec}')):
        _cfg_id = meson_file.strip('meson_')

        print(meson_file.strip('meson_'))
        peram = load_peram(peram_file, Lt, nvec)
        meson_elemental = load_elemental(meson_file, Lt,nvec, mom='mom_0_0_0', disp='disp')
        peram_back = reverse_perambulator_time(peram,nvec)
        print("Contracting pion")
        pion = np.zeros(Lt, dtype=np.cdouble)
        phi_0 = np.einsum("ij,ab->ijab", gamma[5], meson_elemental[0])
        for t in range(0,Lt):
            phi_t = np.einsum("ij,ab->ijab", gamma[5], meson_elemental[t], optimize='optimal')
            tau = peram[t, :, :, :, :]
            tau_ = peram_back[t, :, :, :, :]
            pion[t] = np.einsum("ijab,jkbc,klcd,lida", phi_t, tau, phi_0, tau_, optimize='optimal')
        pion = pion.real
        pion_h5 = h5_utils.create_dset('pion.h5',key='pion',data=pion)


plt.plot(np.arange(0,Lt), pion, '.', label='pion dist')
plt.yscale('log')
plt.legend()
plt.savefig('96.pdf')


# # Determine the path to the directory
# h5_path = os.path.abspath(os.path.join('..','..', 'exolaunch'))

# # Add the directory to sys.path
# if h5_path not in sys.path:
#     sys.path.append(h5_path)
# h5_path

# # Parameters
# Lt = 96
# numvecs_list = [Lt]
# mom_keys = {0: 'mom_-1_0_0', 1: 'mom_-2_0_0', 2: 'mom_-3_0_0', 3: 'mom_0_-1_0', 4: 'mom_0_-2_0', 5: 'mom_0_-3_0', 6: 'mom_0_0_-1', 7: 'mom_0_0_-2', 8: 'mom_0_0_-3', 9: 'mom_0_0_0', 10: 'mom_0_0_1', 11: 'mom_0_0_2', 12: 'mom_0_0_3', 13: 'mom_0_1_0', 14: 'mom_0_2_0', 15: 'mom_0_3_0', 16: 'mom_1_0_0', 17: 'mom_2_0_0', 18: 'mom_3_0_0'}
# mom_keys_inv = {v: k for k, v in mom_keys.items()}
# print(mom_keys_inv)

# disp_keys = {0: 'disp', 1: 'disp_1', 2: 'disp_1_1', 3: 'disp_1_2', 4: 'disp_1_3', 5: 'disp_2', 6: 'disp_2_1', 7: 'disp_2_2', 8: 'disp_2_3', 9: 'disp_3', 10: 'disp_3_1', 11: 'disp_3_2', 12: 'disp_3_3'}
# disp_keys_inv = {v: k for k, v in disp_keys.items()}
# print(disp_keys_inv)
# base_path = Path(h5_path)
# print(base_path)
# configs = list(range(111,231,10))
# print(configs)
# contract.contract(base_path,Lt,32,mom_keys,disp_keys,configs=configs)

