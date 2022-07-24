# Compile models to use in scripts

Compile the model before use by running:

```bash
 python setup.py build_ext --inplace
```

When building a new model edit the setup file 

```python
my_package = Extension('new_model', ['new_model.pyx'],
                     include_dirs=[numpy.get_include()],)
```

and add the package to the setup function: 

```python
setup(ext_modules=cythonize([
                             package1,
                             package2,
                             package3,
                             my_package,
                             ]))
```
