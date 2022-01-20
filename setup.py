
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

# read the version number in as __version__
exec(open('mrpypulse/_version.py').read())

setup(
    name="mrpypulse",
    version=__version__,
    author="Jean-Baptiste Verstraete, Ali Sherzad, Mohammadali Foroozandeh",
    author_email="mohammadali.foroozandeh@chem.ox.ac.uk",
    description=("Mr. PyPulse for Magnetic Resonance pulses"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/foroozandehgroup/mrpypulse",
    packages=find_packages(exclude=["tests"]),
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License"
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=[
        "numpy>=1.17.0",
        "scipy>=1.7.0"
        "matplotlib>=3.3"
    ]
)