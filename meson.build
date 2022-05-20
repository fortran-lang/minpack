project(
  'minpack',
  'fortran',
  version: '2.0.0',
  meson_version: '>=0.55',
  default_options: [
    'default_library=both',
    'buildtype=debugoptimized',
  ],
)
if get_option('api')
  has_cc = add_languages('c', native: false)
elif get_option('python')
  error('Cannot compile Python API, Fortran compiler does not support nested functions')
endif

srcs = files(
  'src/minpack.f90',
)
if get_option('api')
  srcs += files('src/minpack_capi.f90')
endif

minpack_lib = library(
  meson.project_name(),
  sources: srcs,
  install: true,
)

minpack_inc = minpack_lib.private_dir_include()
minpack_dep = declare_dependency(
  link_with: minpack_lib,
  include_directories: [minpack_inc, include_directories('include')],
  variables: {'includedir': meson.current_source_dir() / 'include'},
)

minpack_lic = files(
  'LICENSE.txt',
)

minpack_header = files(
  'include/minpack.h',
)

install_data(
  minpack_lic,
  install_dir: get_option('datadir')/'licenses'/meson.project_name()
)

if get_option('api')
  install_headers(
    minpack_header,
  )
endif

module_id = meson.project_name()
meson.add_install_script(
  find_program(files('config'/'install-mod.py')),
  get_option('includedir') / module_id,
)

pkg = import('pkgconfig')
pkg.generate(
  minpack_lib,
  description: 'Solver for nonlinear equations and least squares problems',
  subdirs: ['', module_id],
)

# add examples
subdir('examples')

# add the testsuite
subdir('test')

if get_option('python')
  subdir('python/minpack')
endif
