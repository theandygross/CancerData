#!/usr/bin/env python

from distutils.core import setup

setup(name='CancerData',
      version='0.1',
      description='Some processing and storage tools for cancer genomics data.',
      author='Andrew Gross',
      author_email='the.andrew.gross@gmail.com',
      url='http://andy-gross.flavors.me',
      package_dir = {'': 'src'},
      packages=['Data', 'Processing'],
     )