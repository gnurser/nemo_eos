.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/nemo_eos.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/nemo_eos
    .. image:: https://readthedocs.org/projects/nemo_eos/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://nemo_eos.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/nemo_eos/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/nemo_eos
    .. image:: https://img.shields.io/pypi/v/nemo_eos.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/nemo_eos/
    .. image:: https://img.shields.io/conda/vn/conda-forge/nemo_eos.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/nemo_eos
    .. image:: https://pepy.tech/badge/nemo_eos/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/nemo_eos
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/nemo_eos
    .. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
       :alt: Project generated with PyScaffold
       :target: https://pyscaffold.org/

********
nemo_eos
********


    f2py'd versions of NEMO Fortran routines used to calculate density, thermal coefficient of expansion etc that can be called by a python script or by ipython or jupyter lab. Currently requires to be installed into an already existing  conda/mamba python environment. If this environment does not exist, setup as described in  Pre-installation, below.


Pre-installation
================
Here we assume that we wish to discard any  existing conda/mamba python setup (normally best). So, firstly remove any such setup

.. code-block :: bash

   cd
   rm -rf miniconda3 anaconda mambaforge
   wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
   sh Miniforge3-Linux-x86_64.sh

Follow prompts; say yes to modifying ``.bashrc`` so that base environment is activated on login

Then create environment that you will use for diagnosis and into which the nemo_eos package will go. I suggest xarray and ipython; you may wish to add others.

.. code-block :: bash

   mamba create -n big python=3.10 xarray ipython

and move into this environment

.. code-block :: bash

   mamba activate big
   


Installation
============

First go into parent directory into which you want to download the package.
  
.. code-block :: bash

   cd /path/to/parent

Then activate big environment if not previously done

.. code-block :: bash

   mamba activate big

Then git clone the source from github and move into new directory

.. code-block :: bash

   git clone https://github.com/JMMP-Group/nemo_eos
   cd nemo_eos

Build and install

.. code-block :: bash

   dev.bash

This installs two shared libraries that calculate the equation of state as used by NEMO. Many of the routines are accelerated by openMP threads.

- ``nemo_rho_omp``        :  old versions of eos in terms of potential tempearture before v4.2, following Jackett and McDougall (1994).
- ``nemo_rho_teos10_omp`` :  new version of eos (both the teos10 and the new implementation of the old expression in terms of potential tempearture.




check works OK and explore
==========================

.. code-block :: bash

   ipython

.. code-block :: python

   Python 3.10.13 | packaged by conda-forge | (main, Dec 23 2023, 15:36:59) [Clang| 16.0.6 ]
   Type 'copyright', 'credits' or 'license' for more information
   IPython 8.22.1 -- An enhanced Interactive Python. Type '?' for help.

   In [1]: from nemo_eos.nemo_rho_omp import eos

   In [2]: eos.grav
   Out[2]: array(9.80665016)

Type in ``eos.<TAB>`` to get list of functions, options, constants etc. Then select routine with cursor keys

.. code-block :: python

   In [3]: eos.rho_mn4
    alpha_beta4        alpha_beta_n4      grav               rho_mn4            s00                sigma_n4          
    alpha_beta8        drho0_mn4          neos               rho_mn8            set_eos            sigma_n8          
    alpha_beta8_nomask drho0_mn8          rau0               rn_alpha           set_eos_threads    theta00           
    alpha_beta_04      fillvalue          rho000             rn_beta            sigma_n                              
    instance

Then when desired routine is selected, type <RET> and other options will disappear

.. code-block :: python

   In [3]: eos.rho_mn4

Then type ?<RET>, and ipython will give description of how to call the routine:

.. code-block :: python

   In [3]: eos.rho_mn4?
   Call signature: eos.rho_mn4(*args, **kwargs)
   Type:           fortran
   String form:    <fortran function rho_mn4>
   Docstring:     
   rho = rho_mn4(fillvalue,mask,theta,s,depth,[n])

   Wrapper for ``rho_mn4``.

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


Check ``eos.rho_mn4``, routine for calculating in-situ density. Check value is ``rho = 1060.93298 kg/m**3`` for ``p=10000 dbar``, ``theta = 40 deg celcius``, ``S=40 psu``

.. code-block :: python

   In [4]: eos.rho_mn4(1.e10, False, 40.0, 40.0, 1.e4)
   Out[4]: array([60.93299], dtype=float32)

Check speed of ``eos.rho_mn4`` for typical number of points for 1 and 4 OpenMP threads

.. code-block :: python

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

   In [8]: timeit rho = eos.rho_mn4(fillvalue,mask,T4,S4,depth4)
   1.05 ms ± 5.82 µs per loop (mean ± std. dev. of 7 runs, 1,000 loops each)

   In [9]: eos.set_eos_threads(4)
   
   In [10]: timeit rho = eos.rho_mn4(fillvalue,mask,T4,S4,depth4)
   320 µs ± 28.6 µs per loop (mean ± std. dev. of 7 runs, 1,000 loops each)   













