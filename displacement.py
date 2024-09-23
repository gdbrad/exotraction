'''representations of derivatives to characterize displacements
there are phases associated with operator displacement
# single pion
#isospinSnk = { (('m',0,2,2),) : 1.0 }
#isospinSrc = { (('m',0,2,2),) : 1.0 }
Let us call a single pion at zero total momentum and zero displacemnt pi_T1
At displacement 1, nabla_pi_T1

'''

pi_src =
pi_snk =



# compute a 2 by 2 two-point correlation matrix
twopt_matrix = twopoint_matrix([op_pi, op_pi2], e, p, list(range(128)), 128)
twopt_matrix = twopt_matrix.real
print(
    backend.arccosh(
        (backend.roll(twopt_matrix[0, 0], -1, 0) + backend.roll(twopt_matrix[0, 0], 1, 0)) / twopt_matrix[0, 0] / 2
    )
)
print(
    backend.arccosh(
        (backend.roll(twopt_matrix[1, 1], -1, 0) + backend.roll(twopt_matrix[1, 1], 1, 0)) / twopt_matrix[1, 1] / 2
    )
)