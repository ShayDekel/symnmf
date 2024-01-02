from setuptools import Extension, setup

module = Extension('symnmfmodule', sources=['symnmf.c', 'symnmfmodule.c'])
setup(
    name='symnmfmodule',
    version='6.9',
    description='This is the greatest project ever made!',
    ext_modules=[module])
