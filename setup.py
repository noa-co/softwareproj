from setuptools import Extension, setup

module = Extension("mykmeanssp",
                   sources=[
                       'spkmeans.c',
                       'spkmeansmodule.c',
                       'sp_io_utils.c',
                       'sp_kmeansalho.c'
                   ])
setup(name='mykmeanssp',
      version='1.0',
      description='Py wrapper for kmeanspp and spkmeans algorithm',
      ext_modules=[module])