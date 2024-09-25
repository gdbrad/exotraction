import h5py
import numpy as np
from typing import List,Dict 
import os
import argparse
import datetime
from ingest_data import load_elemental, load_peram, reverse_perambulator_time
import gamma as gamma
from displacement import displacement_map
from displacement import deriv_names


Insertion = {}
# pion operator with no derivative operators 
Insertion['pi_A1'] = {'gamma':gamma.gamma[5],
                           'deriv':deriv_names.IDEN,
                           'projection': 'A1',
                           'mom':'mom_0_0_0'}
# pion operators with a lattice derivative operator 
# naming: gamma_structre x deriv operator _ O_h irrep 

Insertion['pixnabla_T1']= {'gamma':gamma.gamma[5],
                           'deriv':deriv_names.nabla,
                           'projection': 'T1',
                           'mom':'mom_0_0_0'}
Insertion['pixD_T2'] = {'gamma':gamma.gamma[5] @ gamma.gamma[4],
                           'deriv':deriv_names.D,
                           'projection': 'T2',
                           'mom':'mom_0_0_0'}

Insertion['pixB_T1'] = {'gamma':gamma.gamma[5],
                           'deriv':deriv_names.B,
                           'projection': 'T1',
                           'mom':'mom_0_0_0'}

def displace(disp,u,F,length,mu):
    '''returns sign and coeff'''
    deriv_type = str()
    mu = displacement_map(disp)
    if mu[0]!=mu[1]:
        assert deriv_type=='mixed'


def operator_displacement(name:str,
                          disp:List[int],
                          src:bool,
                          snk:bool,
                          length=1):
    '''O (deriv operator) x gamma construction of a meson operator'''
    ins = Insertion[name]
    if len(disp)==1:
        # length is always = 1
        deriv_operator == ins['deriv']

        displacement(u, F, length, mu) - displacement(u, F, -length, mu)

        # apply first derivative to the right onto src 


        
        
        
     

def contract_displacement(
        cfg_id,
        num_vecs:int,
        num_tsrcs:int,
        peram_dir, 
        meson_dir,
        op_name: List[str],
        disp:List[str],
        Lt:int,
        h5_group, 
        show_plot=False):
    
    peram_file = None
    peram_filename = f"peram_{num_vecs}_cfg{cfg_id}.h5"
    for file in os.listdir(peram_dir):
        if file == peram_filename:
            peram_file = os.path.join(peram_dir, file)
            break

    meson_file = None
    meson_filename = f"meson-{num_vecs}_cfg{cfg_id}.h5"
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
    
    # Load perambulator and meson elemental
    # Once the τ have been computed and stored, the correlation
    # of any source and sink operators can be computed a posteriori.
    nop = len(op_name)
    pion = np.zeros(Lt,nop, dtype=np.cdouble)  # Shape (96, 200) for each tsrc
    peram = load_peram(peram_file, Lt, num_vecs, num_tsrcs)
    peram_back = reverse_perambulator_time(peram)
    for tsrc in range(num_tsrcs):
        for t in range(Lt):
            tau = peram[tsrc, t, :, :, :, :]
            tau_ = peram_back[tsrc, t, :, :, :, :]
            # now loop over displacements
            for _disp in disp:
                meson_elemental = load_elemental(meson_file, Lt, num_vecs, mom='mom_0_0_0', disp=_disp)
                # construct right source with covariant derivative operator eg. nabla,B,D
                _ins = Insertion[op_name]
                phi_0 = np.einsum("ij,ab->ijab", _ins['gamma'], meson_elemental[0])
                for _src in range(nop):
                    for _snk in range(nop):
                        phi_src = phi_0[_src]
                        # apply first right derivative nabla 
                        if len(disp==1):
                            gamma_src = np.einsum("ij,xkj,kl->xil", gamma.gamma[4]@ gamma.gamma[5], phi_src[0].conj(), gamma.gamma[4]@ gamma.gamma[5])
                        phi[0] = operator_displacement(disp=disp)
                        if len(disp==2):

                        phi_t = np.einsum("ij,ab->ijab", gamma.gamma[4]@ gamma.gamma[5], meson_elemental[t], optimize='optimal')

                        contracted_result = np.einsum("ijab,jkbc,klcd,lida", phi_t, tau, phi_0, tau_, optimize='optimal')
                        
                        # Store the contracted result in the pion array (Lt, 200)
                        pion[t] = contracted_result
        
        
        pion = pion.real
        h5_group.create_dataset(f'tsrc_{tsrc}/cfg_{cfg_id}', data=pion)

        if show_plot:
            plt.plot(np.arange(Lt), pion[:, 0], '.', label=f'Pion Distribution (first column) - tsrc {tsrc}, cfg {cfg_id}')
            plt.yscale('log')
            plt.legend()
            plt.savefig(f'pion-{cfg_id}-tsrc-{tsrc}-{num_vecs}-{datetime.datetime.today()}.pdf')

    print(f"Cfg {cfg_id} processed successfully.")
    return True

def main(cfg_ids, num_vecs, num_tsrcs,disp, task_id,show_plot=False):
    h5_path = os.path.abspath('/p/scratch/exotichadrons/exolaunch')
    Lt = 96  
    peram_dir = os.path.join(h5_path, 'perams_sdb', f'numvec{num_vecs}', f'tsrc-{num_tsrcs}')
    meson_dir = os.path.join(h5_path, 'meson_sdb', f'numvec{num_vecs}')



    #return 2pt corr matrix 
    twopt_matrix = contract_displacement(cfg_ids,num_vecs,num_tsrcs,[op_pi, op_pi2], meson_dir, peram_dir, list(range(128)), 128)
    twopt_matrix = twopt_matrix.real

    h5_output_file = f'pion_2pt_nvec_{num_vecs}_tsrc_{num_tsrcs}_task{task_id}.h5'
    with h5py.File(h5_output_file, "w") as h5f:
        h5_group = h5f.create_group('pion_000')
        for cfg_id in cfg_ids:
            try:
                processed = process_configuration(cfg_id, num_vecs, num_tsrcs, peram_dir, meson_dir,disp, Lt, h5_group, show_plot)
                if not processed:
                    print(f"Skipping configuration {cfg_id} file is missing")
            except FileNotFoundError as e:
                print(e)

        print(f"All cfgs processed & saved to {h5_output_file}.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process peram and meson files for multiple configurations.")
    parser.add_argument('--cfg_ids', type=str, required=True, help="List of configuration IDs to process")
    parser.add_argument('--nvec', type=int, required=True, help="Number of eigenvectors")
    parser.add_argument('--ntsrc', type=int, required=True, help="Number of source time slices")
    parser.add_argument('--disp', type=int, required=True, help="displacement")
    parser.add_argument('--plot', action='store_true', help="Show plot of pion distribution")
    parser.add_argument('--task', type=int, required=True, help="SLURM array task ID or unique identifier for this run")


    args = parser.parse_args()
    # if args.cfg_ids is None:
    #     cfg_ids = list(range(11, 1992, 10))  # Generate 11, 21, 31, ..., 1991
    # else:
    cfg_ids =  [int(cfg) for cfg in args.cfg_ids.split(',')]

    main(cfg_ids=cfg_ids, num_vecs=args.nvec, num_tsrcs=args.ntsrc, show_plot=args.plot,disp=args.disp,task_id=args.task)