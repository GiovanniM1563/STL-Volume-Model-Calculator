from setuptools import setup, find_packages

setup(
    name='mypackage',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'math',
        'struct',
        'sys',
        'numpy',
        'trimesh',
        'matplotlib',
        'mpl_toolkits',
    ]
)