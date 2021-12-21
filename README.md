# CAPG - Comprehensive Allopolyploid Genotyper

This software is implemented for genotyping tetraploids' targeted genome regions using screened reads aligned to the whole genome references.
It has primarily been designed as a standalone executable.
Skip to [Installation](#installation) to see how.

# Table of Contents
1. [Prerequisites](#prerequisites)
1. [Installation](#installation)
1. [Tutorial](#tutorial)
1. [Input Files](#input)
1. [Output](#output)
1. [Troubleshooting](#troubleshooting)
1. [Command-Line Options](#options)
1. [How to Cite](#cite)
1. [Acknowledgements](#acknowledgements)
1. [Contact](#contact)

# Prerequisites <a name = "prerequisites" />
- Rmath, the [R Standalone Math Library](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library).  Often, the Rmath library (libRmath.a or libRmath.so for Linux or libRmath.dylib for MacOS) will be installed with R, but not always.  Here are some other locations for the library.
	- r-mathlib on [Ubuntu](https://ubuntu.com/) and [Debian](https://www.debian.org/)
	- libRmath on [Fedora](https://ubuntu.com/), [CentOS](https://centos.org/), [Mageia](https://www.mageia.org/en/), and [Mandriva](https://www.openmandriva.org/)
	- Or if all else fails, you can install the Rmath standalone library from the repository [https://github.com/statslabs/rmath](https://github.com/statslabs/rmath)
# Installation <a name = "installation" />
