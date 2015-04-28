from distutils.core import setup

setup(
    name='pysgrid',
    version='0.0.2',
    author='Andrew Yan',
    author_email='ayan@usgs.gov',
    url='https://github.com/sgrid/pysgrid',
    packages=['pysgrid', 'pysgrid.tests'],
    license='BSD',
    description='Python package for working with staggered gridded data',
    long_description=open('README.md').read(),
    install_requires=['netCDF4 == 1.1.7.1',
                      'numpy == 1.9.2'
                      ],
)