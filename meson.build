project('nemo_eos', 'c',
  version : '0.1',
  license: 'BSD-3',
  meson_version: '>= 1.1.0',
  default_options : ['warning_level=2',
                     'buildtype=debug'],
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


# incdir_f2py = run_command(py,
#     ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
#     check : true
# ).stdout().strip()

#incdir_numpy = meson.get_external_property('numpy-include-dir')#, 'not-given')
incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'

inc_np = include_directories(incdir_numpy, incdir_f2py)

fortranobject_c = incdir_f2py / 'fortranobject.c'

nemo_rho_source = custom_target('nemo_rhomodule.c',
  input : ['src/nemo_eos/fortran/nemo_rho.F90'],
  output : ['nemo_rhomodule.c', 'nemo_rho-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'nemo_rho', '--lower',
             'only:', 'eos_init', 'set_eos_threads', 'get_eos_threads', 'set_eos',
             'get_r0', 'eos_insitu4', 'eos_insitu4_m',
             'eos_sigman4', 'eos_sigman4_m', 'eos_sigma04','eos_sigma04_m',
             'eos_rab4', 'eos_rab4_m', 'eos_rab_ref4', 'eos_rab_ref4_m',
             'eos_pot_from_CT_SA4', 'eos_pen4', 'eos_insitu04', 'eos_insitu04_m',
             ':']
)
py.extension_module('nemo_rho',
  [
    'src/nemo_eos/fortran/nemo_rho.F90', nemo_rho_source
  ],
  fortranobject_c,
  include_directories: inc_np,
  dependencies : py_dep,
  fortran_args :  ['-fcheck=all', '-fopenmp'],
  link_args : ['-lgomp'],
  install : true,
  # install_rpath: '/Users/agn/miniforge3/envs/big/lib',
  subdir: 'nemo_eos'
)
# nemo_rho_teos10_omp_source = custom_target('nemo_rho_teos10_ompmodule.c',
#   input : ['src/nemo_eos/fortran/nemo_rho_teos10_omp.F90'],
#   output : ['nemo_rho_teos10_ompmodule.c', 'nemo_rho_teos10_omp-f2pywrappers2.f90'],
#   command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'nemo_rho_teos10_omp', '--lower']
# )
# py.extension_module('nemo_rho_teos10_omp',
#   [
#     'src/nemo_eos/fortran/nemo_rho_teos10_omp.F90', nemo_rho_teos10_omp_source
#   ],
#   fortranobject_c,
#   include_directories: inc_np,
#   dependencies : py_dep,
#   fortran_args : ['-fopenmp'],
#   link_args : ['-lgomp'],
#   install : true,
#   # install_rpath: '/Users/agn/miniforge3/envs/big/lib',
#   subdir: 'nemo_eos'
#)
python_sources = [
  'src/nemo_eos/__init__.py']
py.install_sources(python_sources, subdir: 'nemo_eos') 


#-mtune=native -march=native
#-Wall -Wextra -Wimplicit-interface -fPIC -O3 -march=native -ffast-math -funroll-loops
