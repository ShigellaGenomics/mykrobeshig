# MykroShig

This package parses Mykrobe predict results for *Shigella sonnei* and *Shigella flexneri*. It is adapted from the [sonneityping tool](https://github.com/katholt/sonneityping) originally developed for *Shigella sonnei*.

Mykrobe v0.9.0+ can identify input genomes and assign them to hierarchical genotypes based on detection of single nucleotide variants (SNVs), and report known mutations in the quinolone-resistance determining region (QRDR) of genes *gyrA* and *parC*.

This package parses the JSON files output by Mykrobe (one per genome) and tabulates the results in a single tab-delimited file.

## Requirements
* Python 3.7+
* No external dependencies (uses only Python standard library)

## Installation

```bash
git clone https://github.com/ShigellaGenomics/mykroshig.git
cd mykroshig
pip install .
```

## Quickstart

```bash
mykroshig --jsons mykrobe_results/*.json --prefix results_mykrobe_parsed
```

## Usage (including running Mykrobe)

### Install Mykrobe
First, install Mykrobe (v0.9.0+) as per the instructions on the [Mykrobe github](https://github.com/Mykrobe-tools/mykrobe).

Once Mykrobe is installed, make sure you run the following two commands to ensure you have the most up-to-date panels for genotyping:
```bash
mykrobe panels update_metadata
mykrobe panels update_species all
```

You can check what version of the scheme is currently loaded in your Mykrobe installation via:
```bash
mykrobe panels describe
```

### Run `mykrobe predict` on each genome

Example command (using Illumina reads):
```bash
mykrobe predict --sample SAMPLE_NAME --species flexneri --format json --out SAMPLE_NAME.json --seq reads_1.fastq.gz reads_2.fastq.gz
```

* For Oxford Nanopore reads, add the flag `--ont` to your command.
* For full details on all Mykrobe options, please see [the Mykrobe documentation](https://github.com/Mykrobe-tools/mykrobe).

### Parse Mykrobe output with MykroShig

Once Mykrobe results are generated, use `mykroshig` to parse them. MykroShig automatically detects whether samples have been typed using the *S. flexneri* or *S. sonnei* schemes from the Mykrobe output.

**Arguments**
* `--jsons`: JSON files output from `mykrobe predict`
* `--prefix`: prefix for output file
* `--force`: If you used `--force` with your Mykrobe call, to enforce a genotype call even when the species and phylogroup calls are not above the internal Mykrobe thresholds, you can force the parser to output this forced genotype for you. Otherwise, genomes with no phylogroup or no species Mykrobe calls will be skipped

```bash
mykroshig --jsons mykrobe_results/*.json --prefix results_mykrobe_parsed
```

## Output
The output table will be named *prefix*_predictResults.tsv, and will be in tab-delimited format.

Explanation of columns in the output:
* **genome**: sample ID
* **species**: species call, one of:
  * *S. flexneri* - if `species coverage` is >85.5
  * *S. sonnei* - if `species coverage` is >90
  * 'Unknown E. coli/Shigella' - if `species coverage` not met, but Mykrobe determines it belongs to `Ecoli_Shigella` genus, then it is likely an E. coli or Shigella genome but not one with a matching genotyping scheme
  * 'Not E. coli/Shigella' - Mykrobe was unable to detect the *uidA* gene, therefore this genome is unlikely to be *E. coli* or *Shigella*
  * 'forced call' - if `--force` provided, then no species is entered here and 'forced call' is used instead
* **final genotype**: final genotype call from Mykrobe
* **name**: human readable alias for genotype. Pulled from the data/alleles_*.txt files. Only used for *S. sonnei*.
* **scheme**: scheme used for typing the genome in Mykrobe (one of sonnei or flexneri)
* **confidence**: measure of confidence in the final genotype call, summarising read support across all levels in the hierarchy (lineage, clade, subclade, etc)
  * _strong_ - high quality calls (quality of '1' reported by Mykrobe) for ALL levels;
  * _moderate_ - reduced confidence for ONE node (Mykrobe quality of '0.5', but still with >50% read support for the marker allele), high quality calls for ALL OTHERS;
  * _weak_ - low quality for one or more nodes (Mykrobe quality of '0.5' and <50% read support OR Mykrobe quality of '0').
* **num QRDR**: Total number of mutations detected in the quinolone-resistance determining regions (QRDR) of genes _gyrA_ and _parC_
* **parC_S80I, gyrA_S83L, gyrA_S83A, gyrA_D87G, gyrA_D87N, gyrA_D87Y**: calls for each individual QRDR mutation. 0 indicates mutation is absent, 1 indicates mutation is present.
* **phylogroup_coverage**: Percent coverage to the *uidA* marker for the E. coli/Shigella grouping, as detected by Mykrobe.
* **species_coverage**: Percent coverage to the species markers for either *S. sonnei* or *S. flexneri*, as detected by Mykrobe. Determines the call in `species` column for the output, based on the thresholds outlined above.
* **lowest support for genotype marker**: For any markers in the final genotype call that do not have a Mykrobe quality of '1', this column reports the percentage of reads supporting the marker allele at the most poorly supported marker
* **poorly supported markers**: Lists any markers in the final genotype call that do not have Mykrobe quality of '1'. Markers are separated by ';', values in brackets represent the quality call from Mykrobe, followed by the read depth at the alternate / reference alleles.
* **max support for additional markers**: For any markers detected that are incongruent with the final genotype call, this column reports the percentage of reads supporting the marker allele at the best supported additional marker.
* **additional markers**: Lists any markers that are incongruent with the final genotype call. Markers are separated by ';', and the format is identical to column _poorly supported markers_.
* **node support**: A list of all markers in the final genotype call with their Mykrobe quality calls (1, 0.5, or 0) and the read depths at the marker allele / reference allele.

### Unexpected results
As recombination patterns may vary in either *S. flexneri* or *S. sonnei*, if you see anything reported in the "additional markers" column, it is worth investigating further whether you have contaminated sequence data (which would produce low-quality calls with low read support at additional markers) or genuine recombination/mixed infections. In such cases it may also be illuminating to investigate the original Mykrobe output file for further information.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

## Citation

If you use the sonnei scheme for typing, please cite: Hawkey, J. et al. Global population structure and genotyping framework for genomic surveillance of the major dysenteric pathogen, Shigella sonnei. Nat Commun 12, 2684 (2021). https://doi.org/10.1038/s41467-021-22700-4

If you use the flexneri scheme for typing, please cite: TBD
