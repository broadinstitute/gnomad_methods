"""Setup script."""

import re

import setuptools


with open("README.md", "r") as readme_file:
    long_description = readme_file.read()

# Hail is deliberately excluded from install_requires: users install the Hail
# version that matches their cluster. Match the name only so that a version
# specifier on the requirements.txt line (e.g. ``hail!=0.2.139``) is still excluded.
install_requires = []
with open("requirements.txt", "r") as requirements_file:
    for req in (line.strip() for line in requirements_file):
        if not req or req.startswith("#"):
            continue
        if re.split(r"[<>=!~\[; ]", req, maxsplit=1)[0] != "hail":
            install_requires.append(req)


setuptools.setup(
    name="gnomad",
    version="0.8.2",
    author="The Genome Aggregation Database",
    author_email="gnomad@broadinstitute.org",
    description="Hail utilities for the Genome Aggregation Database",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/broadinstitute/gnomad_methods",
    packages=setuptools.find_namespace_packages(include=["gnomad.*"]),
    project_urls={
        "Documentation": "https://broadinstitute.github.io/gnomad_methods/",
        "Source Code": "https://github.com/broadinstitute/gnomad_methods",
        "Issues": "https://github.com/broadinstitute/gnomad_methods/issues",
        "Change Log": "https://github.com/broadinstitute/gnomad_methods/releases",
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
    ],
    python_requires=">=3.9",
    install_requires=install_requires,
)
