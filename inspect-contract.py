# Some Sanity checks...
print(np.sum(np.abs(prop[0, :, :, :, :, :] - prop1)))
print(np.sum(np.abs(prop_back[0, :, :, :, :, :] - prop_back1)))
print(np.sum(np.abs(meson[:, mom_keys_inv['mom_0_0_0'], disp_keys_inv['disp'], :, :] - meson1)))
exit()
np.set_printoptions(precision=2, threshold=5)
print(prop[0, 16, 0, 0, :, :])
print()
print(prop_full[16, 0, 1, 3, :, :])
print()
print(prop_back[0, 16, 1, 3, :, :])
print()
print(prop_full[16, 0, 0, 2, :, :] - prop_back[0, 16, 0, 2, :, :])

print(np.sum(np.abs(prop_full[16, 0, :, :, :, :] - prop_back[0, 16, :, :, :, :])))
print(4*4*20*20, 4*4*20*20 * 1e-08)

exit()


((g4 u) g d)' g4 -> d' g' (g4 u)' g4 -> d' g' u' g4' g4

print(np.transpose(np.conj(gamma[4]))*gamma[4])
print(np.transpose(np.conj(gamma[1])) - gamma[1])
exit(1)