project(
  'minpack',
  'c',
  version: '2.0.0',
  meson_version: '>=0.55',
  default_options: [
    'buildtype=debugoptimized',
  ],
)
install = true

minpack_dep = dependency(
  meson.project_name(),
  version: '>=@0@'.format(meson.project_version()),
  fallback: [meson.project_name(), '@0@_dep'.format(meson.project_name())],
  default_options: [
    'default_library=static',
    'api=true',
    'python=false',
  ],
)
minpack_header = files('include'/'_minpack.h')

subdir('minpack')
