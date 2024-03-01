project('nemo_eos', 'c',
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
nemo_eos_source = custom_target('nemo_rho_ompmodule.c',
  input : ['nemo_eos/tempdir/nemo_rho_omp.f90'],  # .f so no F90 wrappers
  output : ['nemo_rho_ompmodule.c', 'nemo_rho_omp-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'nemo_rho_omp', '--lower']#, '-lgomp']
)
inc_np = include_directories(incdir_numpy, incdir_f2py)

py.extension_module('nemo_rho_omp',
  [
    'nemo_eos/tempdir/nemo_rho_omp.f90', nemo_eos_source
  ],
  incdir_f2py / 'fortranobject.c',
  include_directories: inc_np,
  dependencies : py_dep,
  link_args : ['-lgomp'],
  install : true,
  subdir: 'nemo_eos'
)
python_sources = [
  'nemo_eos/__init__.py']
py.install_sources(python_sources, subdir: 'nemo_eos') 



