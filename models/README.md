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
                             #full_model,
                             small_model,
                             small_model_glu,
                             my_package,
                             ]))
```

##  Note on the full model
By default, we don't compile the full model as this requires 
substantial computational resources and might not be suitable
on every system. If you choose to use the "full_model" uncomment 
the full_model in the setup function by deleting the respective # 
in the setup.py file.