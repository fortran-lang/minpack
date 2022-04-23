tests = [
  'chkder',
  'hybrd',
  'hybrj',
  'lmder',
  'lmdif',
  'lmstr',
]

foreach t : tests
  test(
    t,
    executable(
      'test-@0@'.format(t),
      sources: files('test_@0@.f90'.format(t.underscorify())),
      dependencies: minpack_dep,
    ),
    suite: 'unit',
  )
endforeach

if get_option('api')
  test(
    'c-api',
    executable(
      'c-tester',
      sources: files('api/tester.c'),
      dependencies: minpack_dep,
    ),
    suite: 'api',
  )
endif
