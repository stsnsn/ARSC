from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="astur",
    version="0.0.4",
    packages=find_packages(),
    install_requires=["biopython>=1.79"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={
        "console_scripts": [
            "astur=ASTUR.main:main",
        ],
    },
    author="Satoshi Nishino",
    author_email="satoshi-nishino@g.ecc.u-tokyo.ac.jp",
    description="Compute ARSC (N/C/S) from protein fasta files",
    url="https://github.com/stsnsn/ASTUR",
    license="GPL-2.0",
)