# Variant Parser for virtual panels
A parser for selecting a shortlist of variants from an Alamut file.
By specifying a list of candidate genes generated from HPO / Phenomizer, the Alamut input is parsed and a standard ruleset applied in order to reduce the number of variants to a shortlist of potential pathogenic ones. The Alamut input is generated from an WES VCF file.

Created by Garan Jones g.jones@exeter.ac.uk

## Installation

Clone the repository to a suitable directory, initiallizing Git if required;

```bash
	git init
	git clone https://github.com/pasted/variant_parser.git
	bundle install
```
## Usage

Help list of commands

```bash
	ruby variant_parser.rb --help
```

General use, with output to an Excel file
```bash
	ruby variant_parser.rb -v <example_alamut_file.txt> -g <example_gene_list.txt> -p <sample_id> -e

```

```
Options:
  -v, --variants=<s>            Filepath to Alamut file to parse.
  -f, --fake-exome-depth=<s>    Filepath to fake Exome Depth file to parse.
  -g, --genes=<s>               Filepath to text file with list of valid HGVS gene symbols - one symbol per line.
  -o, --ontology-genes=<s>      Filepath to CSV export from HPO web browser.
  -a, --approved-genes          Check gene symbols against HGNC - only process those which are the current symbol
  -e, --excel-output            Export results to an Excel spreadsheet.
  -p, --proband-sample=<s>      Proband sample ID to include in output
  -c, --clinvar                 Include ClinVar pathogenic candidates (CliniVarClinSignifs: Pathogenic; ClinVarReviewStatus: 3 or 4; ClinVarOrigins: Germline).
  -r, --research                Select variants on research criteria.
  -m, --maf-cutoff=<f>          Set a Minor Allele Frequency (MAF) cutoff - integer out of 1 based on ExAC all populations frequencies. (Default: 0.05)
  -s, --result-prefix           Set the prefix for the result file output, for example Family ID. Defaults to the directory name two levels up from the Alamut variant file.
  -l, --all                     Parse all variants without a genelist
  -h, --help                    Show this message
```

##Selection of candidate genes

[See the Wiki section on using HPO / Phenomizer](https://github.com/pasted/variant_parser/wiki/Variant-parser---Examples-of-use)
