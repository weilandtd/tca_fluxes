from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

full_model = Extension('full_model', ['full_model.pyx'],
                     include_dirs=[numpy.get_include()],
                     extra_compile_args=['-O0',])

small_model = Extension('small_model', ['small_model.pyx'],
                     include_dirs=[numpy.get_include()],)

small_model_glu = Extension('small_model_glu', ['small_model_glu.pyx'],
                     include_dirs=[numpy.get_include()],)

setup(ext_modules=cythonize([
                             #full_model,
                             small_model,
                             small_model_glu,
                             ]))
