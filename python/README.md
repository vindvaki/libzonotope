zonotope.py: A python wrapper for libzonotope_c
===============================================

`zonotope.py` wrappes the functionality of `libzonotope_c.so` in a
python module using
[ctypes](http://docs.python.org/library/ctypes.html).


Installation
------------

The module can be installed globally using
[Distutils](http://docs.python.org/distutils/):

```bash
python setup.py install
```

or using [pip](http://www.pip-installer.org/):

```bash
pip install -e .
```

from this directory. The latter also provides an easy way to uninstall
by running

```bash
pip uninstall zonotope
```


Usage
-----

For more information, see [zonotope.py](zonotope.py), or run

```python
help('zonotope')
```

in your python interpreter.
