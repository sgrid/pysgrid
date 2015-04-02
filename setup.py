from distutils.core import setup

setup(
    name='pysgrid',
    version='0.0.1',
    author='Andrew Yan',
    author_email='ayan@usgs.gov',
    packages=['pysgrid', 'pysgrid.tests'],
    license='BSD',
    description='Python package for working with staggered gridded data',
    long_description=open('README.md').read(),
)