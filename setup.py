#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup, find_packages

install_requires = [
    'pyyaml',
    'setuptools',
    'sqlalchemy',
    'zope.sqlalchemy',
    'psycopg2-binary'
]

version = {}
here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'paudm', 'model', 'version.py')) as fp:
    exec(fp.read(), version)

setup(
    name='paudm.model',
    version=version['__version__'],
    packages=find_packages(),
    namespace_packages=['paudm'],
    install_requires=install_requires,
    description="PAUdm model",
    include_package_data=True,
    zip_safe=True,
)
