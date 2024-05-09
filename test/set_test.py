import timeit
import numpy as np
from nemo_eos.nemo_rho import eosbn2

def set_test(num):
     T = 300*np.random.random_sample(num)
     S = 33. + 7*np.random.random_sample(num)
     depth = 4000*np.random.random_sample(num)
     bottom = 4000*np.random.random_sample(num)
     mask = depth>bottom
     T4, S4, depth4 = [x.astype(np.float32) for x in [T,S,depth]]
     fillvalue = np.float32(1.e10)
     return fillvalue,mask,T4,S4,depth4
     

if __name__ =='__main__':
     rho0 = eosbn2.rho0
     print('initialise TEOS-10')
     eosbn2.eos_init(-1)
     print('check r, r0, rho')
     print("r0 should be 4.59763035 kg /m^3 at depth=1000m=1km")
     r0 = eosbn2.get_r0(1.)
     print(f"r0 is {r0:f} kg /m^3 at depth=1000m=1km")
     print("r should be 28.21993233072 kg /m^3 at CT=3°, SA=35.5 g/kg,"
           "depth=1000m=1km")
     rho, = eosbn2.eos_insitu4(3.,35.5,3000.)
     r = rho - eosbn2.get_r0(3.)
     print(f"r is {r:f}  kg /m^3 at CT=3°, SA=35.5 g/kg, depth=1000m=1km")
     print(f"rho is {rho:f}  kg /m^3 at CT=3°, SA=35.5 g/kg, depth=1000m=1km")
     print(f" rho0 is {eosbn2.rho0:f}")
     print("a=alpha*rho0, b=beta*rho0 should be 0.179646281 kg m−3 K−1, 0.765555368 kg m−3 (g/kg)−1"
           " at CT=10°, SA=30.0 g/kg, depth=1000m=1km")
     alpha, beta = eosbn2.eos_rab4(10.,30.,1000.)
     print(f"alpha*rho0, beta*rho0 are {alpha[0]*rho0:f} kg m−3 K−1, {beta[0]*rho0:f} kg m−3 (g/kg)−1"
           " at CT=10°, SA=30.0 g/kg, depth=1000m=1km")

     print('initialise old_EOS-80')
     eosbn2.eos_init(2)
     print("check rho: should be 60.93298 kg/m**3 for p=10000 dbar, t = 40 deg celcius, s=40")
     rho, =  eosbn2.eos_insitu4(40.,40.,10000.)
     print(f"rho is {rho:f}  kg /m^3 at p=10000 dbar, t = 40 deg celcius, s=40")
     print("alpha/beta, beta should be 0.34763 psu/K, 0.72088 x 1.e-3 psu^{-1} "
           " at CT=10°, SA=40.0 g/kg, depth=4000m")
     alpha, beta = eosbn2.eos_rab4(10.,40.,4000.)
     print(f"alpha/beta, beta  are {alpha[0]/beta[0]:f}, {beta[0]*1.e3}x 1.e-3 psu^{-1} "
           " at CT=10°, SA=40.0 g/kg, depth=4000m")
     
     
     #  Check speed of program. Best done with ipython.
     # from nemo_eos.nemo_rho import eosbn2
     # from set_test import set_test
     # fillvalue,mask,T4,S4,depth4 = set_test(100000)
     # eosbn2.eos_init(-1)
     # eosbn2.set_eos_threads(1)
     # timeit rho = eosbn2.eos_insitu4(T4,S4,depth4)
     # timeit rho = eosbn2.eos_insitu4_m(fillvalue, mask,T4,S4,depth4)
     # timeit alpha, beta = eosbn2.eos_rab4(T4,S4,depth4)
     # timeit alpha, beta = eosbn2.eos_rab4_m(fillvalue, mask,T4,S4,depth4)
     # eosbn2.set_eos_threads(4)
     # timeit rho = eosbn2.eos_insitu4(T4,S4,depth4)
     # timeit rho = eosbn2.eos_insitu4_m(fillvalue, mask,T4,S4,depth4)
     # timeit alpha, beta = eosbn2.eos_rab4(T4,S4,depth4)
     # timeit alpha, beta = eosbn2.eos_rab4_m(fillvalue, mask,T4,S4,depth4)

