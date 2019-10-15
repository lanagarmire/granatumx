# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

try:
    long_description = open("README.rst").read()
except IOError:
    long_description = ""

setup(
    name="zgsea",
    version="0.1.0",
    description="A pip package",
    license="MIT",
    author="Xun Zhu",
    packages=find_packages(),
    package_data={'zgsea': ['data/*']},
    install_requires=[],
    long_description=long_description,
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.6",
    ]
)
