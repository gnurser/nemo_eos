import timeit
import numpy as np
from nemo_eos.nemo_rho import eos

def set_test_values(num):
     T = 300*np.random.random_sample(num)
     S = 33. + 7*np.random.random_sample(num)
     depth = 4000*np.random.random_sample(num)
     bottom = 4000*np.random.random_sample(num)
     mask = depth>bottom
     T4, S4, depth4 = [x.astype(np.float32) for x in [T,S,depth]]
     fillvalue = np.float32(1.e10)
     return fillvalue,mask,T4,S4,depth4

if __name__ =='__main__':
     print('Testing nemo_eos against literature test values for rho, alpha, beta')
     eos.set_eos_threads(4)
     print(70*'=','\n')
     rho0 = eos.rho0
     print(f" rho0 is {eos.rho0:f}\n")
     
     print('initialise TEOS-10')
     print(24*'-','\n')
     eos.eos_init(-1)
     print('check r, r0, rho')
     print("r0 should be 4.59763035 kg /m^3 at depth=1000m=1km")
     r0 = eos.get_r0(1.)
     print(f"r0 is {r0:f} kg /m^3 at depth=1000m=1km")
     print()
     print("r should be 28.21993233072 kg /m^3 at CT=3°, SA=35.5 g/kg,"
           " depth=3000m=3km")
     rho, = eos.eos_insitu4(3.,35.5,3000.)
     r = rho - eos.get_r0(3.)
     print(f"r is {r:f}  kg /m^3 at CT=3°, SA=35.5 g/kg, depth=3000m=3km")
     print(f"rho is {rho:f}  kg /m^3 at CT=3°, SA=35.5 g/kg, depth=3000m=3km")
     try:
          from gsw import density
          print("gsw gives \n rho=", density.rho(35.5,3.,3000.)-1000.)
     except:
          print("gsw not installed. To check against gsw library need to do\n"
             "mamba install gsw")
     print()
     print('check alpha, beta')
     print("a=alpha*rho0, b=beta*rho0 should be 0.179646281 kg m−3 K−1, 0.765555368 kg m−3 (g/kg)−1"
           "\n             at CT=10°, SA=30.0 g/kg, depth=1000m=1km")
     alpha, beta = eos.eos_rab4(10.,30.,1000.)
     print(f"alpha*rho0, beta*rho0 are {alpha[0]*rho0:f} kg m−3 K−1, {beta[0]*rho0:f} kg m−3 (g/kg)−1"
           " \n            at CT=10°, SA=30.0 g/kg, depth=1000m=1km")

     
     print('\n initialise old_EOS-80')
     print(24*'-','\n')
     eos.eos_init(2)
     print("check rho: should be 60.93298 kg/m**3 for p=10000 dbar, t = 40 deg celcius, s=40")
     rho, =  eos.eos_insitu4(40.,40.,10000.)
     print(f"rho is {rho:f}  kg /m^3 at p=10000 dbar, t = 40 deg celcius, s=40")
     print()
     print('check alpha, beta')
     print("alpha/beta, beta should be 0.34763 psu/K, 0.72088 x 1.e-3 psu^{-1} "
           "\n            at t=10°, s=40.0 psu, depth=4000m")
     alpha, beta = eos.eos_rab4(10.,40.,4000.)
     print(f"alpha/beta, beta  are {alpha[0]/beta[0]:f}, {beta[0]*1.e3}x 1.e-3 psu^{-1} "
           "\n            at t=10°, s=40.0 psu, depth=4000m")


     print('initialise new polynomial for EOS-80')
     print(24*'-','\n')
     eos.eos_init(0)
     print('check r, r0, rho')
     print("r0 should be 4.59763035 kg /m^3 at depth=1000m=1km")
     r0 = eos.get_r0(1.)
     print(f"r0 is {r0:f} kg /m^3 at depth=1000m=1km")
     print()
     print("r should be 28.35011066567 kg /m^3 at t=3°, S=35.5 psu,"
           "depth=3000m=3km")
     rho, = eos.eos_insitu4(3.,35.5,3000.)
     r = rho - eos.get_r0(3.)
     print(f"r is {r:f}  kg /m^3 t=3°, S=35.5 psu, depth=3000m=3km")
     print(f"rho is {rho:f}  kg /m^3 at t=3°, S=35.5 psu, depth=3000m=3km")
     print("check rho: should be 60.93298 kg/m**3 for p=10000 dbar, t = 40 deg celcius, s=40")
     rho, =  eos.eos_insitu4(40.,40.,10000.)
     print(f"rho is {rho:f}  kg /m^3 at p=10000 dbar, t = 40 deg celcius, s=40")
     print()
     print('check alpha, beta')
     print("alpha/beta, beta should be 0.34763 psu/K, 0.72088 x 1.e-3 psu^{-1} "
          "\n            at t=10°, s=40.0 psu, depth=4000m")
     alpha, beta = eos.eos_rab4(10.,40.,4000.)
     print(f"alpha/beta, beta  are {alpha[0]/beta[0]:f}, {beta[0]*1.e3}x 1.e-3 psu^{-1} "
           "\n            at t=10°, s=40.0 psu, depth=4000m")
