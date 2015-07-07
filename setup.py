from __future__ import with_statement

from setuptools import setup, find_packages

from pysgrid import __version__

reqs = [line.strip() for line in open('requirements.txt')]


def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name                = 'pysgrid',
    version             = __version__,
    description         = 'Python package for working with staggered gridded data',
    author              = 'Andrew Yan',
    author_email        = 'ayan@usgs.gov',
    url                 = 'https://github.com/sgrid/pysgrid',
    packages            = find_packages(),
    license             = 'BSD',
    long_description    = readme(),
    install_requires    = reqs,
    tests_require       = ['mock', 'nose'],
    classifiers         = [
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering',
        ],
    include_package_data = True,
)
