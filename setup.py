#!/usr/bin/env python
from setuptools import setup, find_packages
import sys, os
import versioneer

def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as f:
        return f.read()

setup(name='revmut',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description="REVertant MUTation finder",
      long_description=read("README.rst"),
      classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Topic :: Software Development :: Libraries :: Python Modules'
      ],
      keywords='Python revertant mutation finder HGVS',
      author='Ino de Bruijn',
      author_email='ino@ino.pm',
      url='https://github.com/inodb/revmut',
      license="MIT",
      include_package_data=True,
      packages=find_packages(exclude=['test*']),
      install_requires=['biopython>=1.65'
                        ],
      entry_points={
          'console_scripts': [
              'revmut = revmut.revertant_mutation_checker:main'
          ]
      },
      )
