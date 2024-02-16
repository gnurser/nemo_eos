project('f2py_examples', 'c',
  version : '0.1',
  license: 'BSD-3',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2'],
)

add_languages('fortran')

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

# nemo_eos_source = custom_target('nemo_eosmodule.c',
#   input : ['tempdir/nemo_rho_omp.f90'],  # .f so no F90 wrappers
#   output : ['nemo_eosmodule.c', 'nemo_eos-f2pywrappers.f','nemo_eos-f2pywrappers.f90'],
#   command : [py, '-m', 'numpy.f2py', '@INPUT@', '--f90flags="-fopenmp -march=native -fPIC -shared -O3"', 'lgomp', '-m', 'nemo_eos', '--lower']
# )
nemo_eos_source = custom_target('nemo_eosmodule.c',
  input : ['src/tempdir/nemo_rho_omp.f90'],  # .f so no F90 wrappers
  output : ['nemo_eosmodule.c', 'nemo_eos-f2pywrappers2.f90'],#,'nemo_eos-f2pywrappers.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'nemo_eos', '--lower']#, '-lgomp']
)
inc_np = include_directories(incdir_numpy, incdir_f2py)

py.extension_module('nemo_eos',
  [
    'src/tempdir/nemo_rho_omp.f90', nemo_eos_source
  ],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  link_args : ['-lgomp'],
  install : true
)

