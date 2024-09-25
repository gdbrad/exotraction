import h5py
import numpy as np
import os 
from typing import Iterable,List,Dict,Union
from pathlib import Path
from gamma import gamma

def contract_pion_t(t: int, meson_elemental: np.ndarray, prop: np.ndarray, prop_back: np.ndarray, gamma: dict[int, np.ndarray]) -> np.cdouble:
    phi_0 = np.einsum("ij,ab->ijab", gamma[5], meson_elemental[0])
    phi_t = np.einsum("ij,ab->ijab", gamma[5], meson_elemental[t])
    tau = prop[t, :, :, :, :]
    tau_ = prop_back[t, :, :, :, :]
    p = np.einsum("ijab,jkbc,klcd,lida", phi_t, tau, phi_0, tau_)
    print(t, p)
    return p

def contract_vector_t(t: int, meson_elemental: np.ndarray, prop: np.ndarray, prop_back: np.ndarray, g_idx: int, gamma: dict[int, np.ndarray]) -> np.cdouble:
    phi_0 = np.einsum("ij,ab->ijab", gamma[g_idx], meson_elemental[0])
    phi_t = np.einsum("ij,ab->ijab", np.condjugate(np.transpose(gamma[g_idx])), meson_elemental[t])
    tau = prop[t, :, :, :, :]
    tau_ = prop_back[t, :, :, :, :]
    p = np.einsum("ijab,jkbc,klcd,lida", phi_t, tau, phi_0, tau_)
    print(t, p)
    return p

def contract(
            base_path: Path,
            Lt:int,
            num_vecs:int,
            mom_keys:Dict,
            disp_keys:Dict,
            configs:Union[int, Iterable[int], None] = None)-> Dict[int, Dict[int, np.ndarray]]:
    data = {num_vecs: {}}
    perams_path = os.path.join(base_path,'perams_sdb',f'numvec{num_vecs}')
    elemental_path = os.path.join(base_path,'meson_sdb',f'numvec{num_vecs}')
    output_path = os.path.join(base_path, 'pion_h5', f'numvec{num_vecs}')
    os.makedirs(output_path, exist_ok=True)

    if configs is None:
        configs = range(11, 1992, 10)
    elif isinstance(configs, List):
        configs = configs
    for config in configs:
        peram_file = os.path.join(perams_path, f'peram_{num_vecs}_cfg{config}.h5')
        elemental_file = os.path.join(elemental_path, f'meson-{num_vecs}_cfg{config}.h5')
            
        if os.path.isfile(peram_file) & os.path.isfile(elemental_file):
            print(f"Starting contraction for cfg {config}")
            # continue
        
        with h5py.File(peram_file, 'r') as peram_data:
            peram = np.zeros((Lt, Lt, 4, 4, num_vecs, num_vecs), dtype=np.cdouble)
            for t_source_idx in range(Lt):
                t_source_data = peram_data[f't_source_{t_source_idx}']
                for t_slice_idx in range(Lt):
                    t_slice_data = t_source_data[f't_slice_{t_slice_idx}']
                    for spin_src_idx in range(4):
                        spin_src_data = t_slice_data[f'spin_src_{spin_src_idx}']
                        for spin_snk_idx in range(4):
                            spin_snk_data = spin_src_data[f'spin_sink_{spin_snk_idx}']
                            peram[t_source_idx, t_slice_idx, spin_src_idx, spin_snk_idx, :, :] = \
                                spin_snk_data['real'][:] + spin_snk_data['imag'][:] * 1j

        with h5py.File(elemental_file, 'r') as elemental_data:
            elemental = np.zeros((Lt,  len(mom_keys), len(disp_keys),num_vecs, num_vecs), dtype=np.cdouble)
            for t_slice_idx in range(0, Lt):
                t_slice_data = elemental_data[f't_slice_{t_slice_idx}']
                for mom_idx in range(len(mom_keys.keys())):
                    mom_data = t_slice_data[mom_keys[mom_idx]]
                    for disp_idx in range(0, len(disp_keys)):
                        disp_data = mom_data[disp_keys[disp_idx]]
                        elemental[t_slice_idx, mom_idx, disp_idx, :, :] = \
                    disp_data['real'][:] + disp_data['imag'][:] * 1j
               
                        elemental[t_slice_idx, mom_idx, disp_idx, :, :] = \
                        disp_data['real'][:] + disp_data['imag'][:] * 1j
        mom_keys_inv = {v: k for k, v in mom_keys.items()}
        disp_keys_inv = {v: k for k, v in disp_keys.items()}
    
        meson_elemental = np.squeeze(elemental[:, mom_keys_inv['mom_0_0_0'], disp_keys_inv['disp'], :, :])
        pion = np.zeros((Lt,Lt), dtype=np.cdouble)
        for t_src in range(0,Lt):
            phi_0 = np.einsum("ij,ab->ijab", DP.g5, meson_elemental[t_src])

            for t in range(0,Lt):
                phi_t = np.einsum("ij,ab->ijab", DP.g5, meson_elemental[t])
                tau = np.squeeze(peram[t_src, t, :, :, :, :])
                tau_ = np.squeeze(peram[t, t_src, :, :, :, :])
            
                pion[t_src,t] = np.einsum("ijab,jkbc,klcd,lida", phi_t, tau, phi_0, tau_)
                print(pion[t_src,t].shape, pion[t_src,t])
        pion = pion.real
        output_file = os.path.join(output_path, f'pion_nvec{num_vecs}_cfg{config}.h5')
        with h5py.File(output_file, 'w') as f5:
            f5.create_dataset('pion', data=pion)
        
        print(f"Pion data written to {output_file}")

    return data