"""Setup configuration for mykroshig package"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="mykroshig",
    version="1.0.0",
    author="Jane Hawkey",
    author_email="jane.hawkey@monash.edu",
    description="Parse Mykrobe predict results for Shigella sonnei and Shigella flexneri",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ShigellaGenomics/mykroshig",
    packages=find_packages(),
    package_data={
        "mykroshig": ["data/*.txt"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bioinformatics",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "mykroshig=mykroshig.parser:main",
        ],
    },
)
