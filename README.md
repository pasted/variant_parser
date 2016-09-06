#Variant Parser
A parser for selecting a shortlist of variants from an Alamut file.
By specifying a list of candidate genes generated from HPO / Phenomizer, the Alamut input is parsed and a standard ruleset applied in order to reduce the number of variants to a shortlist of potential pathogenic ones. The Alamut input is generated from an WES VCF file.


##Installation

Clone the repository to a suitable directory, initiallizing Git if required;

```bash
	git init
	git clone https://github.com/pasted/variant_parser.git
	bundle install
```
##Usage

Help list of commands

```bash
	ruby variant_parser.rb --help
```

General use, with output to an Excel file
```bash
	ruby variant_parser.rb -v example_alamut_file.txt -g example_gene_list.txt -e

```

##Selection of candidate genes

[See the Wiki section on using HPO / Phenomizer](https://github.com/pasted/variant_parser/wiki/Variant-parser---Examples-of-use)
