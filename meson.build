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

minpack_lib = library(
  meson.project_name(),
  sources: files(
    'src/minpack.f90',
    'src/minpack_capi.f90',
  ),
  install: true,
)

minpack_inc = minpack_lib.private_dir_include()
minpack_dep = declare_dependency(
  link_with: minpack_lib,
  include_directories: [minpack_inc, include_directories('include')],
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

install_headers(
  minpack_header,
)

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
