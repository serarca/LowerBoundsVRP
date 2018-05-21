from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy


extensions = [
    Extension('cpp_lower_bounds', ['cpp_lower_bounds.pyx', 'lower_bounds.cpp'],
              include_dirs=[numpy.get_include()],
              extra_compile_args=['-std=c++11'],
              language='c++'
              ),
]

setup(
    ext_modules=cythonize(extensions),
    # extra_compile_args=["-w", '-g'],
)