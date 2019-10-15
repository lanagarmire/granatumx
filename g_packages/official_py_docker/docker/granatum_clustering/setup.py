from setuptools import setup, find_packages

import sys, os

VERSION = '0.2.6'

setup(name='granatum_clustering',
      version=VERSION,
      description="granatum_clustering",
      long_description="""""",
      classifiers=[],
      keywords='granatum clustering',
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
