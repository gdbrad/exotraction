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

from typing import List 
def output(deriv_coeff_index):
    c, n = deriv_coeff_index
    assert isinstance(n, int) and n >= 0
    d = ["dx", "dy", "dz"]
    ret = []
    num = -1
    while n >= 0:
        num += 1
        n = n - pow(3, num)
    n = n + pow(3, num)
    for _ in range(num):
        ret.append(d[n % 3])
        n = n // 3
    if c == 1:
        pass
    elif c == -1:
        ret.insert(0, "-")
    else:
        ret.insert(0, str(c))
    return "".join(ret)

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

deriv_names = {
    IDEN:
    NABLA:
    B:
    D: 
}


print(output([1,11]))
print(derivative(3))

def make_operator(
        name: str,
        row: list,
        coeffs: List[float]
):
    '''insert name 
    row of irrep 
    CG coeffs 
    '''
    pass 
# Construct pion interpolators at src and snk 
# this gets passed to the contraction which forms the 2 pt correlator
pi_A1 = Insertion(GammaName.PI, DerivativeName.IDEN, ProjectionName.A1, momentum_dict)
print(pi_A1[0])
op_pi = Operator("pi", [pi_A1[0](0, 0, 0)], [1])

pixnabla_T1 = Insertion(GammaName.PI, DerivativeName.NABLA, ProjectionName.T1, momentum_dict)
print(pixnabla_T1[0])

pixD_T2 = Insertion(GammaName.PI, DerivativeName.D, ProjectionName.T2, momentum_dict)
print(pixD_T2[0])

pixB_T1 = Insertion(GammaName.PI, DerivativeName.T1, ProjectionName.T1, momentum_dict)
print(pixB_T1[0])




op_pi2 = Operator("pi2", [pi_A1[0](0, 0, 0), b1xnabla_A1[0](0, 0, 0)], [3, 1])

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