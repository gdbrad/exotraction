import numpy as np 

sigma: dict[int, np.ndarray] = {}
gamma: dict[int, np.ndarray] = {}
I = np.array([[1, 0],
              [0, 1]])
sigma[1] = np.array([[0, 1],
                    [1, 0]])
sigma[2] = np.array([[0, -1j],
                    [1j, 0]])
sigma[3] = np.array([[1,  0],
                    [0, -1]])

gamma[1] = np.zeros((4, 4), dtype=np.cdouble)
gamma[1][:2, 2:] = 1j * sigma[1]
gamma[1][2:, :2] = -1j * sigma[1]

gamma[2] = np.zeros((4, 4), dtype=np.cdouble)
gamma[2][:2, 2:] = -1j * sigma[2]
gamma[2][2:, :2] = 1j * sigma[2]

gamma[3] = np.zeros((4, 4), dtype=np.cdouble)
gamma[3][:2, 2:] = 1j * sigma[3]
gamma[3][2:, :2] = -1j * sigma[3]

gamma[4] = np.zeros((4, 4), dtype=np.cdouble)
gamma[4][:2, 2:] = I
gamma[4][2:, :2] = I

gamma[5] = gamma[1] @ gamma[2] @ gamma[3] @ gamma[4]