# Installation

Spl-IsoQuant requires a 64-bit Linux system or Mac OS and Python (3.8 and higher) to be pre-installed on it.
Note that Spl-IsoQuant 2.2.x and earlier versions are not well compatible with Python 3.14.
No special hardware is required. 
You will also need

-   [gffutils](https://pythonhosted.org/gffutils/installation.html)
-   [pysam](https://pysam.readthedocs.io/en/latest/index.html)
-   [pyfaidx](https://pypi.org/project/pyfaidx/)
-   [ssw-py](https://pypi.org/project/ssw-py/)
-   [editdistance](https://pypi.org/project/editdistance/)
-   [biopython](https://biopython.org/)
-   [numba](https://numba.pydata.org/)
-   [minimap2](https://github.com/lh3/minimap2)
-   [samtools](http://www.htslib.org/download/)
-   [STAR](https://github.com/alexdobin/STAR) (optional)

## Installing from GitHub
To obtain Spl-IsoQuant you can download repository and install requirements.
Clone Spl-IsoQuant repository and switch to the latest release:
```bash
git clone https://github.com/algbio/spl-IsoQuant.git
cd spl-IsoQuant
git checkout latest
```
Install requirements:
```bash
pip install -r requirements.txt
```
Typically, package installation takes a few minutes.

You also need [samtools](http://www.htslib.org/download/) and [minimap2](https://github.com/lh3/minimap2) to be in the `$PATH` variable.

## Verifying your installation
To verify Spl-IsoQuant installation run
```bash
splisoquant.py --test
```
to run on toy dataset. This should typically take less than 1 minute. 
If the installation is successful, you will find the following information at the end of the log:
```bash
=== Spl-IsoQuant pipeline finished ===
=== TEST PASSED CORRECTLY ===
```
