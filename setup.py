import setuptools
import os
import glob

with open("README.md", "r") as fh:
    long_description = fh.read()


extra_files = []
for (path, directories, filenames) in os.walk('drop/'):
    directories[:] = [d for d in directories if not (d.startswith('.') or d == 'Data')]
    filenames[:] = [f for f in filenames if 
                    not (f.startswith('.') or f.endswith('.Rproj') or f.endswith('.py'))]
    for filename in filenames:
        extra_files.append(os.path.join('..', path, filename))

setuptools.setup(
    name="drop",
    version="0.0.1",
    author="Michaela Müller, Daniela Andrade Salazar",
    author_email="mumichae@in.tum.de",
    description="Detection of RNA Outlier Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://i12g-gagneurweb.informatik.tu-muenchen.de/gitlab/salazar/drop.git",
    packages=setuptools.find_packages(include=["drop", "wBuild", "snakemake"]),
    package_data={'drop': extra_files},
    include_package_data=True
)
