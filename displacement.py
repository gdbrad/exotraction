'''representations of derivatives to characterize displacements
there are phases associated with operator displacement
# single pion
#isospinSnk = { (('m',0,2,2),) : 1.0 }
#isospinSrc = { (('m',0,2,2),) : 1.0 }
Let us call a single pion at zero total momentum and zero displacemnt pi_T1
At displacement 1, nabla_pi_T1

we need to include all combinations of gamma matrices up to three derivatives.

gauge covariant spatial derivatives are combined with a gamma matrix within a fermion bilinear 


'''
from dataclasses import dataclass
@dataclass
class DispGammaMom:
    gamma =           int         ## The gamma matrix for this  displacement
    displacement =    List[int]    ## The displacement path for this gamma
    mom =             Mom_t       ## Array of momenta to generate 

def disp_gamma_mom():



from typing import List 
def derivative(n: int):
    ret = []
    num = -1
    while n >= 0:
        num += 1
        n = n - pow(3, num)
    n = n + pow(3, num)
    for _ in range(num):
        ret.append(n % 3)
        n = n // 3
    return tuple(ret[::-1])

deriv_scheme = {
    "": [
        [[1, 0]],  # 1
    ],
    R"$\nabla$": [
        [[1, 1]],  # dx
        [[1, 2]],  # dy
        [[1, 3]],  # dz
    ],
     R"$\mathbb{B}$": [
        [[1, 11], [-1, 9]],  # dydz-dzdy
        [[1, 6], [-1, 10]],  # dzdx-dxdz
        [[1, 7], [-1, 5]],  # dxdy-dydx
    ],
    # basis vector for T_2
    R"$\mathbb{D}$": [
        [[1, 11], [1, 9]],  # dydz+dzdy
        [[1, 6], [1, 10]],  # dzdx+dxdz
        [[1, 7], [1, 5]],  # dxdy+dydx
    ],}

class deriv_names:
    IDEN =  ""
    NABLA =  R"$\nabla$"
    B = R"$\mathbb{B}$"
    D = R"$\mathbb{D}$"


def displacement_map(disp:str):
    '''displacement directions from right and left '''
    dist_src_snk = List()
    if disp == 'disp':
        dist_src_snk = [0]
    elif disp == 'disp_1':
        dist_src_snk = [1,0]

    elif disp =='disp_1_1':
        dist_src_snk = [1, 1]

    elif disp == 'disp_1_2':
        dist_src_snk = [1, 2]

    elif disp == 'disp_1_3':
        dist_src_snk = [1, 3]

    elif disp == 'disp_2':
        dist_src_snk = [2, 0]

    elif disp == 'disp_2_1':
        dist_src_snk = [2, 1]

    elif disp == 'disp_2_2':
        dist_src_snk = [2, 2]

    elif disp == 'disp_2_3':
        dist_src_snk = [2 ,3]

    elif disp == 'disp_3':
        dist_src_snk = [3, 0]

    elif disp == 'disp_3_1':
        dist_src_snk = [3, 1]

    elif disp == 'disp_3_2':
        dist_src_snk = [3, 2]

    elif disp == 'disp_3_3':
        dist_src_snk = [3, 3]



print(output([1,11]))
print(derivative(3))

# Construct pion interpolators at src and snk 
# this gets passed to the contraction which forms the 2 pt correlator

# pi_src =
# pi_snk =



# # compute a 2 by 2 two-point correlation matrix
# twopt_matrix = twopoint_matrix([op_pi, op_pi2], e, p, list(range(128)), 128)
# twopt_matrix = twopt_matrix.real
# print(
#     backend.arccosh(
#         (backend.roll(twopt_matrix[0, 0], -1, 0) + backend.roll(twopt_matrix[0, 0], 1, 0)) / twopt_matrix[0, 0] / 2
#     )
# )
# print(
#     backend.arccosh(
#         (backend.roll(twopt_matrix[1, 1], -1, 0) + backend.roll(twopt_matrix[1, 1], 1, 0)) / twopt_matrix[1, 1] / 2
#     )
# )