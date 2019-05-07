# wBuild

Automatic build tool for R Reports

## Documentation

Full documentation is available at https://wbuild.readthedocs.io

## Features

* Supports reproducible research
  * Append R-markdown scripts to the Snakemake pipeline
* Render the R scripts to a structured web page
* Snakemake rules written directly in the header of scripts
  * Dependencies updated automatically!

## Installation

Install using pip:

- `pip install wBuild`

See "Installation" tab in the documentation for more details.

## Getting started
  
* Navigate to an empty directory
* Run `wbuild demo`. This will create a wBuild demo project with various examples
* Explore the files in `Scripts/`
* Run `snakemake` to build the project
* Open `Output/html/index.html` in your web browser
  * There you will find useful and understandable examples of features and operations with wBuild

## Usage

* Navigate to the root of your project (new or existing)
* Run `wbuild init`
* Run `snakemake`

## GitHub

wBuild is an open-source software. The source-code is available at https://github.com/gagneurlab/wBuild.

## Credits

Leonhard Wachutka
