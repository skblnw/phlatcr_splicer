#!/usr/bin/env python3
"""
Setup script for pHLA-TCR Complex Structure Analyzer
"""

from setuptools import setup, find_packages
import os

# Read README for long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read requirements
with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="phlatcr-splicer",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A Python tool for analyzing pHLA-TCR complex structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/phlatcr_splicer",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest",
            "black",
            "flake8",
            "mypy",
        ],
        "docs": [
            "sphinx",
            "sphinx-rtd-theme",
        ],
    },
    entry_points={
        "console_scripts": [
            "phlatcr-analyze=phlatcr_splicer.analyzer:main",
            "mhc-ii-analyze=phlatcr_splicer.mhc_ii_analyzer:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)