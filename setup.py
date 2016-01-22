from __future__ import (absolute_import, division, print_function,
                        with_statement)

import os
from setuptools import setup, find_packages


def extract_version(module='pysgrid'):
    version = None
    fdir = os.path.dirname(__file__)
    fnme = os.path.join(fdir, module, '__init__.py')
    with open(fnme) as fd:
        for line in fd:
            if (line.startswith('__version__')):
                _, version = line.split('=')
                # Remove quotation characters.
                version = version.strip()[1:-1]
                break
    return version


reqs = [line.strip() for line in open('requirements.txt')]


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pysgrid',
      version=extract_version(),
      description='Python package for working with staggered gridded data',
      author='Andrew Yan',
      author_email='ayan@usgs.gov',
      url='https://github.com/sgrid/pysgrid',
      packages=find_packages(),
      license='BSD',
      long_description=readme(),
      install_requires=reqs,
      tests_require=['mock', 'nose'],
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          ],
      include_package_data=True,
      )
