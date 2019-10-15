from setuptools import setup, find_packages

import sys, os

VERSION = '0.0.4'

setup(name='granatum_deeplearning',
      version=VERSION,
      description="granatum_deeplearning",
      long_description="""""",
      classifiers=[],
      keywords='granatum deeplearning',
      author='o_poirion',
      author_email='o.poirion@gmail.com',
      url='',
      license='MIT',
      packages=find_packages(exclude=['examples',
                                      'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[],
      )
