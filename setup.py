from setuptools import setup, find_packages

setup(
    name='VaspTools',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'lxml',
        'pandas',
        'polars',
        'numpy',
        'matplotlib',
        'pymatgen'
    ]
)
