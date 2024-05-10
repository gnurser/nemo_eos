<!-- .. These are examples of badges you might want to add to your README: -->
<!--    please update the URLs accordingly -->

<!--     .. image:: https://api.cirrus-ci.com/github/<USER>/nemo_eos.svg?branch=main -->
<!--         :alt: Built Status -->
<!--         :target: https://cirrus-ci.com/github/<USER>/nemo_eos -->
<!--     .. image:: https://readthedocs.org/projects/nemo_eos/badge/?version=latest -->
<!--         :alt: ReadTheDocs -->
<!--         :target: https://nemo_eos.readthedocs.io/en/stable/ -->
<!--     .. image:: https://img.shields.io/coveralls/github/<USER>/nemo_eos/main.svg -->
<!--         :alt: Coveralls -->
<!--         :target: https://coveralls.io/r/<USER>/nemo_eos -->
<!--     .. image:: https://img.shields.io/pypi/v/nemo_eos.svg -->
<!--         :alt: PyPI-Server -->
<!--         :target: https://pypi.org/project/nemo_eos/ -->
<!--     .. image:: https://img.shields.io/conda/vn/conda-forge/nemo_eos.svg -->
<!--         :alt: Conda-Forge -->
<!--         :target: https://anaconda.org/conda-forge/nemo_eos -->
<!--     .. image:: https://pepy.tech/badge/nemo_eos/month -->
<!--         :alt: Monthly Downloads -->
<!--         :target: https://pepy.tech/project/nemo_eos -->
<!--     .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter -->
<!--         :alt: Twitter -->
<!--         :target: https://twitter.com/nemo_eos -->
<!--     .. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold -->
<!--        :alt: Project generated with PyScaffold -->
<!--        :target: https://pyscaffold.org/ -->

# nemo_eos

f2py'd versions of NEMO Fortran routines used to calculate density,
thermal coefficient of expansion etc that can be called by a python
script or by ipython or jupyter lab. Depending on how it is
initialised, uses : 
(`neos=2`) polynomial TEOS-10 EOS (Roquet et al., 2015),  
(`neos=-1`) MacDougall and Jackett 1994; Macdougall 1987  
 (`neos=0`)  polynomial approximation to EOS-80 (Roquet, pers. comm.)  
 (`neos=1`)  Simplified EOS (vallis, 2006)  
 
 Currently requires to be installed into an already existing  conda/mamba python environment. If this environment does not exist, setup as described in  Pre-installation, below.


## Pre-installation

```bash
	cd  
	rm -rf miniconda3 anaconda mambaforge  
	wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh  
	sh Miniforge3-Linux-x86_64.sh  
```

Follow prompts; say **yes** to modifying `.bashrc` so that base environment is activated on login

Then create environment that you will use for diagnosis and into which the nemo_eos package will go. I suggest xarray and ipython; you may wish to add others.
```
   mamba create -n big python=3.10 xarray ipython
```

and move into this environment

```
   mamba activate big
```


## Installation

First go into parent directory into which you want to download the package.
  
```
    cd /path/to/parent
```

Then activate `big` environment if not previously done

```
    mamba activate big
```

Then git clone the source from github and move into new directory

```
   git clone https://github.com/gnurser/nemo_eos.git
   cd nemo_eos
```

Build and install

```
   dev.bash
```

This installs a shared library `nemo_rho`  that contains the
routines that calculate the equation of state as used by NEMO. Many of
the routines are accelerated by openMP threads.

## Available routines

- `eos_init(neos)`: Set EOS type
- `get_r0` :  For TEOS10 EOS calculate depth-varying offset `r0=rho-r`
  (Roquet et al., 2015)
- `eos_insitu4` : calculate in situ density for 1D arrays of `T, S, depth`
- `eos_insitu4_m`:  calculate in situ density for 1D arrays of `T, S,
  depth`; needs 1D array `mask` and `fill_value`
-  `eos_sigman4` : calculate density referenced to depth `depth_km` km
                    for 1D arrays `T, S`
- `eos_sigman4_m`: masked version of `eos_sigman4`

- `eos_rab_ref4`

Public variables include:

-  Default Boussinesq ocean density `rho0`, EOS type `neos`
-  Coefficients for Vallis's simplified EOS:
	 `rn_lambda1, rn_nu, rn_lambda2, rn_mu1, rn_mu2, rn_a0, rn_b0`
 
            eos_insitu04    eos_pen4        eos_rab_ref4_m                    set_eos        
            eos_insitu04_m  eos_rab4        eos_sigma04     get_eos_threads                    set_eos_threads
                 eos_rab4_m      eos_sigma04_m 

## Test
```
	cd test
	python set_test.py
```

## Explore

```
   ipython

   Python 3.10.13 | packaged by conda-forge | (main, Dec 23 2023, 15:36:59) [Clang| 16.0.6 ]
   Type 'copyright', 'credits' or 'license' for more information
   IPython 8.22.1 -- An enhanced Interactive Python. Type '?' for help.

   In [1]: from nemo_eos.nemo_rho import eos

   In [2]: eos.rho0
   Out[2]: array(1026.)

```

Type in `eos.<TAB>` to get list of functions, options, constants etc. Then select routine with cursor keys

```

   In [3]: eos.eos_insitu4_m
            eos_init        eos_insitu4_m   eos_rab_ref4    eos_sigman4     neos            rn_lambda1      rn_nu          
            eos_insitu04    eos_pen4        eos_rab_ref4_m  eos_sigman4_m   rho0            rn_lambda2      set_eos        
            eos_insitu04_m  eos_rab4        eos_sigma04     get_eos_threads rn_a0           rn_mu1          set_eos_threads
            eos_insitu4     eos_rab4_m      eos_sigma04_m   get_r0          rn_b0           rn_mu2 
```
Then when desired routine is selected, type <RET> and other options will disappear
```

   In [3]: eos.eos_insitu4_m

```
Then type ?<RET>, and ipython will give description of how to call the routine:

```

   In [3]: eos.eos_insitu4_m?
   Call signature: eos.rho_mn4(*args, **kwargs)
   Type:           fortran
   String form:    <fortran function eos_insitu4_m>
   Docstring:     
   rho = eos_insitu4_m(fillvalue,mask,theta,s,depth,[n])

   Wrapper for `eos_insitu4_m`.

   Parameters
   ----------
   fillvalue : input float
   mask : input rank-1 array('b') with bounds (n)
   theta : input rank-1 array('f') with bounds (n)
   s : input rank-1 array('f') with bounds (n)
   depth : input rank-1 array('f') with bounds (n)

   Other Parameters
   ----------------
   n : input int, optional
       Default: shape(mask, 0)

   Returns
   -------
   rho : rank-1 array('f') with bounds (n)

```
Check `eos.eos_insitu4_m`, routine for calculating in-situ density. Check value is `rho = 1060.93298 kg/m**3` for `p=10000 dbar`, `theta = 40 deg celcius`, `S=40 psu`

```

   In [4]: eos.eos_insitu4_m(1.e10, False, 40.0, 40.0, 1.e4)
   Out[4]: array([60.93299], dtype=float32)
```
Check speed of `eos.eos_insitu4_m` for typical number of points for 1 and 4 OpenMP threads

```

   In [4]: import numpy as np
   In [6]: def setup(num):
       ...:     T = 300*np.random.random_sample(num)
       ...:     S = 33. + 7*np.random.random_sample(num)
       ...:     depth = 4000*np.random.random_sample(num)
       ...:     bottom = 4000*np.random.random_sample(num)
       ...:     mask = depth>bottom
       ...:     T4, S4, depth4 = [x.astype(np.float32) for x in [T,S,depth]]
       ...:     fillvalue = np.float32(1.e10)
       ...:     return fillvalue,mask,T4,S4,depth4
       ...:

   In [6]: fillvalue,mask,T4,S4,depth4 = setup(100000)

   In [7]: eos.set_eos_threads(1)

   In [8]: timeit rho = eos.eos_insitu4_m(fillvalue,mask,T4,S4,depth4)
   1.05 ms ± 5.82 µs per loop (mean ± std. dev. of 7 runs, 1,000 loops each)

   In [9]: eos.set_eos_threads(4)
   
   In [10]: timeit rho = eos.eos_insitu4_m(fillvalue,mask,T4,S4,depth4)
   320 µs ± 28.6 µs per loop (mean ± std. dev. of 7 runs, 1,000 loops each)   
```












