import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="drop",
    version="0.0.1",
    author="Michaela Müller, Daniela Andrade Salazar",
    author_email="mumichae@in.tum.de",
    description="Detection of RNA Outlier Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://i12g-gagneurweb.informatik.tu-muenchen.de/gitlab/salazar/drop.git",
    packages=setuptools.find_packages(include=["drop", "wBuild", "snakemake"])
)
