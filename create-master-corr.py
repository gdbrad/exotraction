import sys
import os
import numpy as np
import h5py
import argparse
import matplotlib.pyplot as plt
from contract_gb import contract_pion_t, contract_vector_t
from ingest_data_gb import load_elemental, load_peram, reverse_perambulator_time
import gamma_gb as gamma
import datetime

def process_configuration(cfg_id, num_vecs, tsrc, peram_dir, meson_dir, Lt, h5_group, show_plot=False):
    peram_file = None
    peram_filename = f"peram_32_cfg{cfg_id}.h5"
    for file in os.listdir(peram_dir):
        if file == peram_filename:
            peram_file = os.path.join(peram_dir, file)
            break

    meson_file = None
    meson_filename = f"meson-32_cfg{cfg_id}.h5"
    for file in os.listdir(meson_dir):
        if file == meson_filename:
            meson_file = os.path.join(meson_dir, file)
            break

    if not peram_file:
        print(f"Peram file for configuration {cfg_id} not found in {peram_dir}. Skipping.")
        return False
    if not meson_file:
        print(f"Meson file for configuration {cfg_id} not found in {meson_dir}. Skipping.")
        return False

    print(f"Reading propagator file: {peram_file}")
    print(f"Reading meson elementals file: {meson_file}")
    
    peram = load_peram(peram_file, Lt, num_vecs,num_tsrcs=24)
    meson_elemental = load_elemental(meson_file, Lt, num_vecs, mom='mom_0_0_0', disp='disp')
    peram_back = reverse_perambulator_time(peram)

    # Contract pion
    
    pion = np.zeros(Lt, dtype=np.cdouble)
    phi_0 = np.einsum("ij,ab->ijab", gamma.gamma[5], meson_elemental[0])
    for t in range(Lt):
        phi_t = np.einsum("ij,ab->ijab", gamma.gamma[5], meson_elemental[t], optimize='optimal')
        tau = peram[tsrc,t, :, :, :, :]
        tau_ = peram_back[tsrc,t, :, :, :, :]
        pion[t] = np.einsum("ijab,jkbc,klcd,lida", phi_t, tau, phi_0, tau_, optimize='optimal')

    pion = pion.real
    pion_h5_name = f'pion_{cfg_id}.h5'
    with h5py.File(pion_h5_name, "w") as f:
       f.create_dataset(str(cfg_id), data=pion)


    if show_plot:
        plt.plot(np.arange(Lt), pion, '.', label=f'Pion Distribution - {cfg_id}')
        plt.yscale('log')
        plt.legend()
        plt.savefig(f'pion-{cfg_id}-{num_vecs}-{datetime.datetime.today()}.pdf')

    h5_group.create_dataset(f'cfg_{cfg_id}', data=pion)
    print(f"Configuration {cfg_id} processed successfully.")
    return True

def main(num_configs, num_vecs, tsrc, show_plot=False):
    h5_path = os.path.abspath('/home/grant/jureca_exolaunch')
    Lt = 96  # Lattice temporal extent

    peram_dir = os.path.join(h5_path, 'perams_sdb', f'numvec{num_vecs}', f'tsrc-{tsrc}')
    meson_dir = os.path.join(h5_path, 'meson_sdb', f'numvec{num_vecs}')

    # HDF5 file to store the data for all configurations
    h5_output_file = f'pion_contractions_numvecs_{num_vecs}_tsrc_{tsrc}.h5'
    with h5py.File(h5_output_file, "w") as h5f:
        for tsrc_idx in range(tsrc):
            group_name = f'tsrc_{tsrc_idx}'
            h5f.create_dataset(group_name, data=data[tsrc_idx], compression="gzip")

        h5_group = h5f.create_group('pion_contractions')

        config_ids = range(11, 11 + num_configs * 10, 10)  # Config IDs like 11, 21, 31...

        for cfg_id in config_ids:
            try:
                processed = process_configuration(cfg_id, num_vecs, tsrc, peram_dir, meson_dir, Lt, h5_group, show_plot)
                if not processed:
                    print(f"Skipping configuration {cfg_id} due to missing files.")
            except FileNotFoundError as e:
                print(e)

        print(f"All configurations processed and saved to {h5_output_file}.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process peram and meson files for multiple configurations.")
    parser.add_argument('--num_configs', type=int, required=True, help="Number of configurations to process")
    parser.add_argument('--num_vecs', type=int, required=True, help="Number of eigenvectors")
    parser.add_argument('--tsrc', type=int, required=True, help="Source time slice")
    parser.add_argument('--show_plot', action='store_true', help="Show plot of pion distribution")

    args = parser.parse_args()

    main(num_configs=args.num_configs, num_vecs=args.num_vecs, tsrc=args.tsrc, show_plot=args.show_plot)
