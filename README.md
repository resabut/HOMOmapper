# HOMOmapper: A Tool for Finding Homologous Regions Between Individuals PLINK SNP Data
Author: Joan Escrivà Font  
Version: 1.0
## Description
HOMOmapper is a powerful tool that enables the identification of homologous regions between one individual and a list of references. This homologous regionsand termed ROH (runs of homology).

By using a reference individual, it locates homologous regions in the other individual using a sliding window approach, with the window size determined by the minimum ROH length.
The tool allows users to set the maximum number of mismatches within the window, and generates a table with essential information including the start and end position of ROHs, the chromosome, number of mismatches, and ROH length.
It is applicable for various scenarios, where the reference individual could be a close relative or a population representative.

## Usage

```bash
python3 HOMOmapper.py [-h] -i indiv_ped -r ref_ped [-o output_txt] [-s] [-sl short_len] -l min_roh_length -m max_mismatch [-b]
```
#### Positional Arguments:
- indiv_ped: Individual ped file
- ref_ped: Reference ped file
#### Required Arguments:
- -l min_roh_length, --min-roh-length min_roh_length: Minimum ROH length. It is used for the length of the sliding window (default: 100)
- -m max_mismatch, --max-mismatch max_mismatch: Maximum number of mismatches within the sliding window (default: 5)
#### Optional Arguments:
- -h, --help: Show this help message and exit
- -o output_txt, --out output_txt: Output file (default: out_roh_table.txt)
- -s, --short: Shorten ped files. Useful for testing
- -sl short_len, --short_len short_len: Number of columns (Phenotype + genotype data) to be kept in the shortened ped files. Useful for testing (default: 500)
- -b, --bash: Use bash to sort ped files. Faster but requires UNIX, gawk


## Installation
The main script is a python 3 scrip that requires at least python 3.6.  
The script requires the following packages to be installed:  
    - pandas 1.5.3  
    - tqdm 4.62.3

In UNIX based systems, the following commands can be used to install the required packages:
```bash
pip install pandas==1.5.3
pip install tqdm==4.62.3
```
For the quicker bash sorting, gawk is required. If not installed, it can be installed with the following command:
```bash
sudo apt-get install gawk
```

## Example runs
### Example 1
This dataset contains data for 4 individuals (3 reference individuals and 1 individual to be compared) that are
closely related, they belong to the same population (Tuscany).
It can be used to test the tool and to understand the output as following:
```bash
python3 HOMOmapper.py -i Data/Example1/indiv.ped -r Data/Example1/ref.ped -o Results/Example1/out_roh_table.txt -l 40 -m 5 -b
```
This runs the tool with the following parameters:
- Individual ped file: Data/Example1/indiv.ped
- Reference ped file: Data/Example1/ref.ped
- Output file: Results/Example1/out_roh_table.txt
- Minimum ROH length: 40 (sliding window size)
- Maximum number of mismatches: 5
- Bash sorting: True

### Example 2
This dataset contains data for 5 individuals (4 reference individuals and 1 individual to be compared) that are not
closely related, they belong to different populations (Tuscany, Palestine, Yemen, Hungary and Crete).
It can be used to test the tool and to understand the output as following:
```bash
python3 HOMOmapper.py -i Data/Example2/indiv.ped -r Data/Example2/ref2.ped -o Results/Example2/out_roh_table.txt -l 40 -m 5 -b
```
This runs the tool with the following parameters:
- Individual ped file: Data/Example2/indiv.ped
- Reference ped file: Data/Example2/ref2.ped
- Output file: Results/Example2/out_roh_table.txt
- Minimum ROH length: 40 (sliding window size)
- Maximum number of mismatches: 5
- Bash sorting: True


# Credits
This tool was developed by Joan Escrivà Font.


# Contact
For any questions, please contact Joan Escrivà Font at esfontjo@gmail.com.