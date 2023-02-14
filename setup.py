#!/bin/env python

#######################################################################
#  Copyright (C) 2020 Vinh Tran
#
#  fdog is the python package for feature-aware directed ortholog
#  search. fdog is a free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  fdog is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with fdog.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

from setuptools import setup, find_packages

with open("README.md", "r") as input:
    long_description = input.read()

setup(
    name="fdog",
    version="0.1.8",
    python_requires='>=3.7.0',
    description="Feature-aware Directed OrtholoG search tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Vinh Tran",
    author_email="tran@bio.uni-frankfurt.de",
    url="https://github.com/BIONF/fDOG",
    packages=find_packages(),
    package_data={'': ['*']},
    install_requires=[
        'biopython',
        'tqdm',
        'ete3',
        'six',
        'PyYAML',
        'pyhmmer',
        'pysam',
        'greedyFAS>=1.11.2'
    ],
    entry_points={
        'console_scripts': ["fdog.run = fdog.runSingle:main",
                            "fdogs.run = fdog.runMulti:main",
                            "fdog.setup = fdog.setupfDog:main",
                            "fdog.checkData = fdog.checkData:main",
                            "fdog.addTaxon = fdog.addTaxon:main",
                            "fdog.addTaxa = fdog.addTaxa:main",
                            "fdog.showTaxa = fdog.showTaxa:main",
                            "fdog.setPaths = fdog.setPaths:main",
                            "fdog.mergeOutput = fdog.mergeOutput:main",
                            "fdog.uninstall = fdog.removefDog:main",
                            "fdog.assembly = fdog.fDOGassembly:main",
                            "fdog.mergeAssembly = fdog.mergeAssemblyOutput:main"],
    },
    license="GPL-3.0",
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)
