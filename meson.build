project('nemo_eos', 'c',
  version : '0.1',
  license: 'BSD-3',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2', 'buildtype=release'],
)

add_languages('fortran')
ff = meson.get_compiler('fortran')
# if ff.has_argument('-Wno-conversion')
#   add_project_arguments('-Wno-conversion', language: 'fortran')
# endif

# add_project_arguments('-march':'native', language: 'fortran')
# add_project_arguments('-fopenmp', language: 'fortran')

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

# nemo_eos_source = custom_target('nemo_eosmodule.c',
#   input : ['tempdir/nemo_rho_omp.f90'],  # .f so no F90 wrappers
#   output : ['nemo_eosmodule.c', 'nemo_eos-f2pywrappers.f','nemo_eos-f2pywrappers.f90'],
#   command : [py, '-m', 'numpy.f2py', '@INPUT@', '--f90flags="-fopenmp -march=native -fPIC -shared -O3"', 'lgomp', '-m', 'nemo_eos', '--lower']
# )
nemo_rho_omp_source = custom_target('nemo_rho_ompmodule.c',
  input : ['nemo_eos/tempdir/nemo_rho_omp.f90'],  # .f so no F90 wrappers
  output : ['nemo_rho_ompmodule.c', 'nemo_rho_omp-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'nemo_rho_omp', '--lower']#, '-lgomp']
)
py.extension_module('nemo_rho_omp',
  [
    'nemo_eos/tempdir/nemo_rho_omp.f90', nemo_rho_omp_source
  ],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  link_args : ['-lgomp'],
  install : true,
  install_rpath: '/Users/agn/miniforge3/envs/big/lib',
  subdir: 'nemo_eos'
)
nemo_rho_teos10_omp_source = custom_target('nemo_rho_teos10_ompmodule.c',
  input : ['nemo_eos/tempdir/nemo_rho_teos10_omp.f90'],  # .f so no F90 wrappers
  output : ['nemo_rho_teos10_ompmodule.c', 'nemo_rho_teos10_omp-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'nemo_rho_teos10_omp', '--lower']#, '-lgomp']
)
py.extension_module('nemo_rho_teos10_omp',
  [
    'nemo_eos/tempdir/nemo_rho_teos10_omp.f90', nemo_rho_teos10_omp_source
  ],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  fortran_args : ['-fopenmp'],
  link_args : ['-lgomp'],
  install : true,
  install_rpath: '/Users/agn/miniforge3/envs/big/lib',
  subdir: 'nemo_eos'
)
python_sources = [
  'nemo_eos/__init__.py']
py.install_sources(python_sources, subdir: 'nemo_eos') 


#-mtune=native -march=native
#-Wall -Wextra -Wimplicit-interface -fPIC -O3 -march=native -ffast-math -funroll-loops
