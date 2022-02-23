tests = [
  'hybrd',
  'hybrd1',
  'lmder1',
  'lmdif1',
  'primes',
]

foreach t : tests
  test(
    t,
    executable(
      'example-@0@'.format(t),
      sources: files('example_@0@.f90'.format(t.underscorify())),
      dependencies: minpack_dep,
    ),
    suite: 'example',
  )
endforeach
