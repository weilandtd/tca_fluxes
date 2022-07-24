from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

package1 = Extension('full_model', ['full_model.pyx'],
                     include_dirs=[numpy.get_include()],
                     extra_compile_args=['-O0',])

package2 = Extension('small_model', ['small_model.pyx'],
                     include_dirs=[numpy.get_include()],)

package3 = Extension('small_model_glu', ['small_model_glu.pyx'],
                     include_dirs=[numpy.get_include()],)

setup(ext_modules=cythonize([
                             package1,
                             package2,
                             package3,
                             ]))
