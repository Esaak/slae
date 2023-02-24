import scipy as sp
import numpy as np



def create_test():
    #o = np.eye(3, dtype=float)
    q = np.arange(1, 10).reshape(3,3)

    print(np.array(np.linalg.qr(q)[0]) )
create_test()