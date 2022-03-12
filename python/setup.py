from setuptools import setup

setup(
    cffi_modules=["ffi-builder.py:ffibuilder"],
    package_data={"minpack": ["_libminpack*.so"]},
)

