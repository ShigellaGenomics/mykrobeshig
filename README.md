# mykroshig

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

## Usage

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

### Run mykrobe predict on each genome (Illumina reads, nanopore reads, etc)

Example command (using Illumina reads):
```bash
mykrobe predict --sample SAMPLE_NAME --species flexneri --format json --out SAMPLE_NAME.json --seq reads_1.fastq.gz reads_2.fastq.gz
```

* For Oxford Nanopore reads, add the flag `--ont` to your command.
* For full details on all Mykrobe options, please see [the Mykrobe documentation](https://github.com/Mykrobe-tools/mykrobe).

### Parse Mykrobe output

Once Mykrobe results are generated, use `mykroshig` to parse them:

**Input**
* JSON files output from `mykrobe predict` (`--jsons`)
* Species to parse: either `sonnei` or `flexneri` (`--species`)

**Output**
* tab-delimited file with one row per genome detailing genotype and QRDR mutations (`--prefix`)

**Example commands**

```bash
mykroshig --jsons mykrobe_results/*.json --species flexneri --prefix results_mykrobe_parsed
```


## Example output
The output table will be named *prefix*_predictResults.tsv, and will be in tab-delimited format:

| genome     | species   | final genotype | name      | confidence        | num QRDR | parC_S80I | gyrA_S83L | gyrA_S83A | gyrA_D87G | gyrA_D87N | gyrA_D87Y | lowest support for genotype marker | poorly supported markers      | max support for additional markers | additional markers         | node support                                                                                                                    |
|------------|-----------|----------------|-----------|-------------------|----------|-----------|-----------|-----------|-----------|-----------|-----------|------------------------------------|-------------------------------|------------------------------------|----------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| sampleA | S.flexneri | 1.1.1        | Example_Type      | strong            | 2        | 1         | 1         | 0         | 0         | 0         | 0         |                                    |                               |                                    |                            | lineage1 (1; 97/0); lineage1.1 (1; 120/0); lineage1.1.1 (1; 91/0)                                   |

Explanation of columns in the output:
* **genome**: sample ID
* **species**: species call from Mykrobe (*S. flexneri* or unknown)
* **final genotype**: final genotype call from Mykrobe, using the S. flexneri genotyping scheme
* **name**: human readable alias for genotype, where available
* **confidence**: measure of confidence in the final genotype call, summarising read support across all levels in the hierarchy (lineage, clade, subclade, etc)
  * _strong_ - high quality calls (quality of '1' reported by Mykrobe) for ALL levels;
  * _moderate_ - reduced confidence for ONE node (Mykrobe quality of '0.5', but still with >50% read support for the marker allele), high quality calls for ALL OTHERS;
  * _weak_ - low quality for one or more nodes (Mykrobe quality of '0.5' and <50% read support OR Mykrobe quality of '0').
* **num QRDR**: Total number of mutations detected in the quinolone-resistance determining regions (QRDR) of genes _gyrA_ and _parC_
* **parC_S80I, gyrA_S83L, gyrA_S83A, gyrA_D87G, gyrA_D87N, gyrA_D87Y**: calls for each individual QRDR mutation. 0 indicates mutation is absent, 1 indicates mutation is present.
* **lowest support for genotype marker**: For any markers in the final genotype call that do not have a Mykrobe quality of '1', this column reports the percentage of reads supporting the marker allele at the most poorly supported marker
* **poorly supported markers**: Lists any markers in the final genotype call that do not have Mykrobe quality of '1'. Markers are separated by ';', values in brackets represent the quality call from Mykrobe, followed by the read depth at the alternate / reference alleles.
* **max support for additional markers**: For any markers detected that are incongruent with the final genotype call, this column reports the percentage of reads supporting the marker allele at the best supported additional marker.
* **additional markers**: Lists any markers that are incongruent with the final genotype call. Markers are separated by ';', and the format is identical to column _poorly supported markers_.
* **node support**: A list of all markers in the final genotype call with their Mykrobe quality calls (1, 0.5, or 0) and the read depths at the marker allele / reference allele.

### Unexpected results
As recombination patterns may vary in *S. flexneri*, if you see anything reported in the "additional markers" column, it is worth investigating further whether you have contaminated sequence data (which would produce low-quality calls with low read support at additional markers) or genuine recombination/mixed infections. In such cases it may also be illuminating to investigate the Mykrobe output file for further information.

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

## Citation

If you use the sonnei scheme for typing, please cite: Hawkey, J. et al. Global population structure and genotyping framework for genomic surveillance of the major dysenteric pathogen, Shigella sonnei. Nat Commun 12, 2684 (2021). https://doi.org/10.1038/s41467-021-22700-4

If you use the flexneri scheme for typing, please cite: TBD
