project(
  'minpack',
  'c',
  meson_version: '>=0.53',
  default_options: [
    'buildtype=debugoptimized',
  ],
)
install = true

minpack_dep = dependency('minpack', version: '>=2.0.0')
minpack_header = files('../include/minpack.h')

subdir('minpack')
