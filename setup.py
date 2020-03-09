import setuptools


with open("README.md", "r") as readme_file:
    long_description = readme_file.read()

with open("requirements.txt", "r") as requirements_file:
    install_requires = [line.strip() for line in requirements_file if line.strip()]


setuptools.setup(
    name="gnomad",
    version="0.1.0",
    author="The Genome Aggregation Database",
    author_email="exomeconsortium@gmail.com",
    description="Hail utilities for the Genome Aggregation Database",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/broadinstitute/gnomad_methods",
    packages=setuptools.find_namespace_packages(include=["gnomad.*"]),
    project_urls={
        "Documentation": "https://broadinstitute.github.io/gnomad_methods/",
        "Source Code": "https://github.com/broadinstitute/gnomad_methods",
        "Issues": "https://github.com/broadinstitute/gnomad_methods/issues",
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
    ],
    python_requires=">=3.6",
    install_requires=install_requires,
)
