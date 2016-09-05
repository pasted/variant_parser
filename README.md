#Variant Parser
A parser for selecting a shortlist of variants from an Alamut file.


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
