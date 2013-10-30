#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from setuptools import setup, find_packages

install_requires = [
    'astropysics',
    'matplotlib',
    'numpy',
    'pyfits==3.0.11',
    'pyyaml',
    'scipy',
    'setuptools',
    'sqlalchemy',
    'zope.sqlalchemy',
]

# Read README and CHANGES files for the long description
here = os.path.abspath(os.path.dirname(__file__))
#README  = open(os.path.join(here, 'README.txt')).read()

# Read the version information
execfile(os.path.join(here, 'paudm', 'tools', 'release.py'))

setup(
    name = 'paudm.tools',
    version = __version__, # @UndefinedVariable
    packages = find_packages(),
    namespace_packages = ['paudm'],
    
    install_requires = install_requires,
    
    description = "PAUdm tools",
    #long_description = README,
    #author = 'Pau Tallada Crespí',
    #author_email = 'pau.tallada@gmail.com',
    #maintainer = 'Pau Tallada Crespí',
    #maintainer_email = 'pau.tallada@gmail.com',
    
    include_package_data=True,
    zip_safe=True,
)