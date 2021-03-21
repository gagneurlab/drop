import setuptools
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()
    requirements = [x for x in requirements if not x.startswith("#") and x != ""]

extra_files = []
for (path, directories, filenames) in os.walk('drop/'):
    directories[:] = [d for d in directories if not d.startswith('.')]
    filenames[:] = [f for f in filenames if not f.startswith('.') and not f.endswith('.Rproj')]
    for filename in filenames:
        extra_files.append(os.path.join('..', path, filename))

setuptools.setup(
    name="drop",
    version="1.0.3",
    author="Michaela Müller, Daniela Klaproth-Andrade, Vicente Yépez, Christian Mertes",
    author_email="mumichae@in.tum.de",
    description="Detection of RNA Outlier Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gagneurlab/drop",
    packages=setuptools.find_packages(include=["drop", "drop.*"]),
    entry_points={'console_scripts': ['drop=drop.cli:main']},
    package_data={'drop': extra_files},
    include_package_data=True,
    install_requires=requirements,
)

