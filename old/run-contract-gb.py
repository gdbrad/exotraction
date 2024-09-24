# #         ############ BOOTSTRAP #########################
# #         '''at this point we have a h5 file with dimensions
# #         N_t x t_srcs x N_cfg x Num(nvecs)
# #         96 x 96 x 200 x 4
# #         '''
# #         import bs_utils

# #         ########### AVERAGE ###########################

import sys
import os
import numpy as np 
from pathlib import Path
import argparse
import h5py 
import matplotlib.pyplot as plt
from contract_gb import contract_pion_t, contract_vector_t
from ingest_data_gb import load_elemental, load_peram, reverse_perambulator_time
import h5_utils
import gamma_gb as gamma

def get_elemental(meson_file, Lt, num_vecs, mom:str, disp:str):
    meson_elemental = load_elemental(meson_file, Lt, num_vecs, mom=mom, disp=disp)
    return meson_elemental


def main(cfg_id, num_vecs, tsrc, show_plot=False):
    h5_path = os.path.abspath('/home/grant/jureca_exolaunch')

    Lt = 96  # Lattice temporal extent

    # Define paths
    peram_dir = os.path.join(h5_path, 'perams_sdb', f'numvec{num_vecs}', f'tsrc-{tsrc}')
    meson_dir = os.path.join(h5_path, 'meson_sdb', f'numvec{num_vecs}')


    # Locate peram file for the given cfg_id
    peram_file = None
    for file in os.listdir(peram_dir):
        if file.endswith(f"{cfg_id}.h5"):
            peram_file = os.path.join(peram_dir, file)
            break

    if not peram_file:
        raise FileNotFoundError(f"Peram file for configuration {cfg_id} not found in {peram_dir}")

    # Locate meson file for the given cfg_id
    meson_file = None
    for file in os.listdir(meson_dir):
        if file.endswith(f"{cfg_id}.h5"):
            meson_file = os.path.join(meson_dir, file)
            break

    if not meson_file:
        raise FileNotFoundError(f"Meson file for configuration {cfg_id} not found in {meson_dir}")

    # Load peram and meson data
    peram = load_peram(peram_file, Lt, num_vecs)
    peram_back = reverse_perambulator_time(peram)

    
    meson_elemental = get_elemental(meson_file, Lt, num_vecs, mom='mom_0_0_0', disp='disp')

    # Contract pion
    pion = np.zeros(Lt, dtype=np.cdouble)
    phi_0 = np.einsum("ij,ab->ijab", gamma.gamma[5], meson_elemental[0])
    for t in range(Lt):
        phi_t = np.einsum("ij,ab->ijab", gamma.gamma[5], meson_elemental[t], optimize='optimal')
        tau = peram[t, :, :, :, :]
        tau_ = peram_back[t, :, :, :, :]
        pion[t] = np.einsum("ijab,jkbc,klcd,lida", phi_t, tau, phi_0, tau_, optimize='optimal')

    pion = pion.real
    import datetime
    today = datetime.datetime.today()

    if show_plot:
        plt.plot(np.arange(Lt), pion, '.', label='Pion Distribution')
        plt.yscale('log')
        plt.legend()
        plt.savefig(f'pion-{num_vecs}-{today}.pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process a single configuration for peram and meson files.")
    parser.add_argument('--cfg_id', type=int, required=True, help="Configuration ID to process")
    parser.add_argument('--num_vecs', type=int, required=True, help="Number of eigenvectors")
    parser.add_argument('--tsrc', type=int, required=True, help="Source time slice")
    parser.add_argument('--show_plot', action='store_true', help="Show plot of pion distribution")

    args = parser.parse_args()

    main(cfg_id=args.cfg_id, num_vecs=args.num_vecs, tsrc=args.tsrc, show_plot=args.show_plot)
