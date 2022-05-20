cc = meson.get_compiler('c')
ext_module = '_libminpack'

pymod = import('python')
python = pymod.find_installation(
  get_option('python_version'),
  modules: [
    'cffi',
  ],
)
python_dep = python.dependency(required: true)

# Python's CFFI is horrible in working with preprocessor statements,
# therefore, we have to preprocess the header before passing it to the ffibuilder
minpack_pp = configure_file(
  command: [
    cc,
    '-I@0@'.format(
      minpack_dep.get_variable(
        pkgconfig: 'includedir',
        internal: 'includedir',
      ).split().get(0)
    ),
    '-DMINPACK_CFFI=1',
    '-E',
    '@INPUT@',
  ],
  input: minpack_header,
  output: '@0@.h'.format(ext_module),
  capture: true,
)

# This is the actual out-of-line API processing of the ffibuilder
minpack_cffi_srcs = configure_file(
  command: [python, files('..'/'ffi-builder.py'), '@INPUT@', '@BASENAME@'],
  input: minpack_pp,
  output: '@BASENAME@.c',
)

# Actual generation of the Python extension
minpack_pyext = python.extension_module(
  ext_module,
  minpack_cffi_srcs,
  dependencies: [minpack_dep, python_dep],
  install: true,
  subdir: 'minpack',
)

python.install_sources(
  '__init__.py',
  'exception.py',
  'library.py',
  'typing.py',
  'test_library.py',
  subdir: 'minpack',
)
