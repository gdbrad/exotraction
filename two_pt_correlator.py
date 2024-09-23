from typing import Literal,NamedTuple,List,Tuple,Iterable
from sympy import Add, Mul, Dummy, Symbol, S, simplify, sqrt
import opt_einsum
from dataclasses import dataclass
import numpy as np
import h5py as h5
import pathlib

import gamma
DP = gamma.DiracPauli

#  this should ingest data from input file 
def C_2pt(
    operators: List[Operator],
    elemental: np.ndarray,
    perambulator: np.ndarray,
    Lt: int,
    timeslices: Iterable[int],
    num_vecs:List[int],
    sum_source_t: bool,
    meson_file: pathlib.Path,
    peram_file: pathlib.Path,
    num_nabla: int,
    mom_list: List[Tuple[int]],
    disp_list: List[Tuple[int]]
    )->np.ndarray:
    '''
    C(t, 0) = Tr[phi(t)tau(t,0)phi(0)tau(0,t)]
    Calculate the two-point correlation function for a given set of operators.

    Parameters:
    - operators: List of Operator objects defining the interpolating fields.
    - elemental: Meson elemental object .
    - perambulator: perambulator (quark propagator) data.
    - timeslices: Iterable of time slices at which the correlation function is to be evaluated.
    - Lt: Temporal extent of the lattice.
    - numvecs: Number of eigenvectors used in the calculation (will be <<< distillation basis)
    - is_sum_over_source_t: return all t source correlator at axis = 0, default is True.

    Returns:
    - A NumPy array of shape (Nop, Lt) containing the two-point correlation function
      values for each operator and timeslice.

    '''
    # Nop = len(operators)
    Nt = len(timeslices)

    # result = np.zeros(Nop,Nt,Lt, dtype=np.cdouble)
    result = np.zeros(Nt,Lt, dtype=np.cdouble)

    peram,meson = load_data(peram_file,meson_file,Lt,num_vecs,mom_list,disp_list)
    phis = collect_elementals(operators,elemental,num_vecs)
    for mom in mom_list:
        for disp in disp_list:
            elemental = np.squeeze(meson[:, mom_list[mom], disp_list[disp], :, :])
            for it,t in enumerate(Nt):
                tau = np.squeeze(peram[0, t, :, :, :, :])
                #tmp = peram[t, :, :, :, :num_vecs, :num_vecs] 
                # if we generated perams for nvec= num eigvecs, we could do this 
                # but we are computing perams for different values of nvecs <<< eigs_numvecs 
                tmp = np.squeeze(peram[t, 0, :, :, :, :])
                tau_back = opt_einsum.contract("ii,tjiba,jj->tijab", DP.g_5, tmp.conj(), DP.g_5)
                for op_idx in range(Nop):
                    phi = phis[op_idx]
                    result[op_idx, it] = opt_einsum.contract(
                        "tijab,xjk,xtbc,tklcd,yli,yad->t",
                        tau_back,
                        phi[0],
                        np.roll(phi[t], -t, 1),
                        tau,
                        phi[0],
                        np.conjugate(phi[t][:, t]),
            )
                    phi_0 = np.einsum("ij,ab->ijab", gamma_and_derivatives.g_5, elemental[0])
                    phi_t = np.einsum("ij,ab->ijab", gamma_and_derivatives.g_5, elemental[t])

                # for a in range(4):
                #     for b in range(4):
                #         tau_[a,b,:,:] = np.transpose(np.conjugate(tau[a,b,:,:]))
                result[t] = np.einsum("ijab,jkbc,klcd,lida", phi_t, tau, phi_0, tau_)
                print(result[t].shape, result[t])
    pion = pion.real