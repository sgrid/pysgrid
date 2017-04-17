from __future__ import (absolute_import, division, print_function,
                        with_statement)

import os
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def run_tests(self):
        # Import here, cause outside the eggs aren't loaded.
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)


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


install_requires = [line.strip() for line in open('requirements.txt')]
tests_require = [line.strip() for line in open('requirements-dev.txt')]


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
      install_requires=install_requires,
      tests_require=tests_require,
      cmdclass=dict(test=PyTest),
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Developers',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: BSD License',
                   'Operating System :: POSIX :: Linux',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering',
                   ],
      include_package_data=True,
      zip_safe=False
      )
