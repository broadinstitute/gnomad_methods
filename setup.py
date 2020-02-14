import setuptools


with open("README.md", "r") as readme_file:
    long_description = readme_file.read()

with open("requirements.txt", "r") as requirements_file:
    install_requires = [line.strip() for line in requirements_file if line.strip()]


setuptools.setup(
    name="gnomad_hail",
    version="0.0.1",
    author="The Genome Aggregation Database",
    author_email="exomeconsortium@gmail.com",
    description="Hail utilities for the Genome Aggregation Database",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/macarthur-lab/gnomad_hail",
    packages=setuptools.find_packages(),
    project_urls={
        "Documentation": "https://macarthur-lab.github.io/gnomad_hail/",
        "Source Code": "https://github.com/macarthur-lab/gnomad_hail",
        "Issues": "https://github.com/macarthur-lab/gnomad_hail/issues",
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.6",
    install_requires=install_requires,
)
