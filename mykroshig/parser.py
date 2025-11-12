#     This script parses either sonnei or flexneri genotyping results output by Mykrobe predict (JSON format)
#     Adapted from sonneityping script by Kat Holt and Jane Hawkey
#
#     Copyright (C) 2025  Jane Hawkey
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

import json
import sys
import csv
from argparse import ArgumentParser

# Species-specific configuration
SPECIES_CONFIG = {
    'flexneri': {
        'species_name': 'Shigella_flexneri',
        'display_name': 'S. flexneri',
        'qrdr_mutations': ['parC_S80I', 'gyrA_S83L', 'gyrA_S83A', 'gyrA_D87G', 'gyrA_D87N', 'gyrA_D87Y'],
        'coverage_threshold': 80,
        'phylo_group': 'Ecoli_Shigella',
    },
    'sonnei': {
        'species_name': 'Shigella_sonnei',
        'display_name': 'S. sonnei',
        'qrdr_mutations': ['parC_S80I', 'parC_S80R', 'gyrA_S83L', 'gyrA_S83A', 'gyrA_D87G', 'gyrA_D87N', 'gyrA_D87Y'],
        'coverage_threshold': 80,
        'phylo_group': 'Ecoli_Shigella',
    }
}


def get_arguments():
    parser = ArgumentParser(description='Parse mykrobe predict JSON files for Shigella sonnei or Shigella flexneri')

    parser.add_argument('--jsons', required=True, nargs='+', help='JSON files output from mykrobe predict')
    parser.add_argument('--prefix', required=True, help='prefix for output files')
    parser.add_argument('--species', required=True, choices=['sonnei', 'flexneri'],
                        help='Species to parse - either sonnei or flexneri')
    parser.add_argument('--alleles', required=False, help='Optional alleles file for sonnei genotyping')

    return parser.parse_args()


def extract_qrdr_info(genome_data, genome_name, spp_call, config):
    """Extract QRDR mutation information for the genome"""
    qrdr_out_dict = {}
    qrdr_possible = config['qrdr_mutations']

    # can only extract qrdr info if sample matches the species, otherwise set 'NA' for all calls
    if spp_call == config['species_name']:
        try:
            qrdr_data = genome_data["susceptibility"]["ciprofloxacin"]
        except KeyError:
            qrdr_out_dict['genome'] = genome_name
            for allele in qrdr_possible:
                qrdr_out_dict[allele] = 'NA'
            qrdr_out_dict['num QRDR'] = 'NA'
            return qrdr_out_dict
    else:
        qrdr_out_dict['genome'] = genome_name
        for allele in qrdr_possible:
            qrdr_out_dict[allele] = 'NA'
        qrdr_out_dict['num QRDR'] = 'NA'
        return qrdr_out_dict

    # set up the list of calls in this genome
    qrdr_calls = []
    # only need to parse if the predict is R - if it's S all values will be 0
    if qrdr_data["predict"] == "R":
        calls = qrdr_data["called_by"]
        # loop through each mutation
        for mutation in calls:
            # get the mutation name and add it to our list of calls
            qrdr_calls.append(mutation.split('-')[0])
    
    # now create our output dictionary for this genome
    qrdr_out_dict['genome'] = genome_name
    num_qrdr_calls = 0
    for allele in qrdr_possible:
        if allele in qrdr_calls:
            qrdr_out_dict[allele] = 1
            num_qrdr_calls += 1
        else:
            qrdr_out_dict[allele] = 0
    # add column with total number of qrdr calls
    qrdr_out_dict['num QRDR'] = num_qrdr_calls
    return qrdr_out_dict


def inspect_calls(full_lineage_data):
    """Inspect genotype calls and determine the best genotype with confidence scores"""
    genotype_details = full_lineage_data['calls_summary']
    genotype_list = list(genotype_details.keys())
    
    best_score = 0
    best_genotype = None
    
    # loop through each genotype
    for genotype in genotype_list:
        max_score = genotype_details[genotype]['tree_depth']
        actual_score = sum(list(genotype_details[genotype]['genotypes'].values()))
        
        if actual_score == max_score:
            best_score = actual_score
            best_genotype = genotype
        elif actual_score < max_score and actual_score > best_score:
            best_score = actual_score
            best_genotype = genotype

    best_calls = genotype_details[best_genotype]['genotypes']
    best_calls_vals = list(best_calls.values())
    poorly_supported_markers = []
    lowest_within_genotype_percents = {}
    final_markers = []
    
    for level in best_calls.keys():
        call_details = full_lineage_data['calls'][best_genotype][level]
        
        if call_details:
            call_details = call_details[list(call_details.keys())[0]]
            ref = call_details['info']['coverage']['reference']['median_depth']
            alt = call_details['info']['coverage']['alternate']['median_depth']
            
            try:
                percent_support = alt / (alt + ref)
            except ZeroDivisionError:
                percent_support = 0
            
            marker_string = level + ' (' + str(best_calls[level]) + '; ' + str(alt) + '/' + str(ref) + ')'
            final_markers.append(marker_string)
        else:
            lowest_within_genotype_percents[0] = level
            marker_string = level + ' (0)'
            poorly_supported_markers.append(marker_string)
            final_markers.append(marker_string)

        if best_calls[level] < 1:
            lowest_within_genotype_percents[percent_support] = level
            poorly_supported_markers.append(marker_string)

    # Determine confidence level
    if best_calls_vals.count(0) == 0 and best_calls_vals.count(0.5) == 0:
        confidence = 'strong'
        lowest_support_val = ''
    elif best_calls_vals.count(0) == 0 and best_calls_vals.count(0.5) == 1:
        if min(lowest_within_genotype_percents.keys()) > 0.5:
            confidence = 'moderate'
        else:
            confidence = 'weak'
        lowest_support_val = round(min(lowest_within_genotype_percents.keys()), 3)
    else:
        confidence = 'weak'
        lowest_support_val = round(min(lowest_within_genotype_percents.keys()), 3)

    # Handle additional incongruent markers
    non_matching_markers = []
    non_matching_supports = []
    
    if len(genotype_list) > 1:
        genotype_list.remove(best_genotype)
        for genotype in genotype_list:
            other_calls = genotype_details[genotype]['genotypes']
            for call in other_calls.keys():
                if call not in best_calls.keys():
                    call_info = full_lineage_data['calls'][genotype][call]
                    if call_info:
                        call_info = call_info[list(call_info.keys())[0]]
                        ref_depth = call_info['info']['coverage']['reference']['median_depth']
                        alt_depth = call_info['info']['coverage']['alternate']['median_depth']
                        
                        if alt_depth >= 1:
                            percent_support = alt_depth / (alt_depth + ref_depth)
                            non_matching_supports.append(percent_support)
                            marker_string = call + ' (' + str(other_calls[call]) + '; ' + str(alt_depth) + '/' + str(ref_depth) + ')'
                            non_matching_markers.append(marker_string)

    if len(non_matching_supports) > 0:
        max_non_matching = round(max(non_matching_supports), 3)
    else:
        max_non_matching = ''

    return best_genotype, confidence, lowest_support_val, poorly_supported_markers, max_non_matching, non_matching_markers, final_markers


def extract_lineage_info(lineage_data, genome_name, config):
    """Extract lineage and species information"""
    spp_data = lineage_data['species']
    spp_call = list(spp_data.keys())[0]
    phylo_grp_data = lineage_data['phylo_group']
    phylo_grp_call = list(phylo_grp_data.keys())[0]
    
    # if spp is unknown, then this is not the target species
    if spp_call == "Unknown":
        if phylo_grp_call != 'Unknown':
            phylo_grp_percentage = phylo_grp_data[phylo_grp_call]['percent_coverage']
        else:
            phylo_grp_percentage = 'NA'
        out_dict = {
            'genome': genome_name, 
            'species': f'not {config["display_name"]}', 
            'final genotype': 'NA',
            'confidence': 'NA',
            'phylogroup_coverage': phylo_grp_percentage,
            'species_coverage': 'NA',
            'lowest support for genotype marker': '', 
            'poorly supported markers': '', 
            'node support': '', 
            'max support for additional markers': '', 
            'additional markers': ''
        }
        return out_dict, spp_call
    else:
        # if it is the target species, then get the percentage
        spp_percentage = spp_data[config['species_name']]["percent_coverage"]
        phylo_grp_percentage = phylo_grp_data[config['phylo_group']]['percent_coverage']
        
        # if the percentage is below threshold, then exit this function
        if spp_percentage < config['coverage_threshold']:
            out_dict = {
                'genome': genome_name, 
                'species': f'not {config["display_name"]}', 
                'final genotype': 'NA',
                'confidence': 'NA',
                'phylogroup_coverage': phylo_grp_percentage,
                'species_coverage': spp_percentage,
                'lowest support for genotype marker': '', 
                'poorly supported markers': '', 
                'node support': '', 
                'max support for additional markers': '', 
                'additional markers': ''
            }
            return out_dict, "Unknown"

    # Continue with target species
    spp_percentage = spp_data[config['species_name']]["percent_coverage"]
    phylo_grp_percentage = phylo_grp_data[config['phylo_group']]['percent_coverage']
    lineage_out_dict = {'genome': genome_name}

    try:
        genotype_calls = lineage_data['lineage']['lineage']
    except KeyError:
        out_dict = {
            'genome': genome_name, 
            'species': config['display_name'], 
            'final genotype': 'uncalled',
            'confidence': 'NA',
            'phylogroup_coverage': phylo_grp_percentage,
            'species_coverage': spp_percentage,
            'lowest support for genotype marker': '', 
            'poorly supported markers': '', 
            'max support for additional markers': '', 
            'additional markers': '', 
            'node support': ''
        }
        return out_dict, spp_call
    
    # if there are no calls, populate with none
    if len(genotype_calls) == 0:
        lineage_out_dict['final genotype'] = 'uncalled'
        lineage_out_dict['confidence'] = 'NA'
        lineage_out_dict['lowest support for genotype marker'] = ''
        lineage_out_dict['poorly supported markers'] = ''
        lineage_out_dict['max support for additional markers'] = ''
        lineage_out_dict['additional markers'] = ''
        lineage_out_dict['node support'] = ''
    else:
        # inspect calls for genotype(s)
        best_genotype, confidence, lowest_support_val, poorly_supported_markers, non_matching_support, non_matching_markers, final_markers = inspect_calls(lineage_data['lineage'])
        lineage_out_dict['final genotype'] = best_genotype
        lineage_out_dict['confidence'] = confidence
        lineage_out_dict['lowest support for genotype marker'] = lowest_support_val
        lineage_out_dict['poorly supported markers'] = '; '.join(poorly_supported_markers)
        lineage_out_dict['max support for additional markers'] = non_matching_support
        lineage_out_dict['additional markers'] = '; '.join(non_matching_markers)
        lineage_out_dict['node support'] = '; '.join(final_markers)
    
    # add species info
    lineage_out_dict['species'] = config['display_name']
    lineage_out_dict['phylogroup_coverage'] = phylo_grp_percentage
    lineage_out_dict['species_coverage'] = spp_percentage

    return lineage_out_dict, spp_call


def main():
    args = get_arguments()
    config = SPECIES_CONFIG[args.species]
    
    results_data = []

    # read in json files
    for json_file in args.jsons:
        try:
            with open(json_file) as f:
                myk_result = json.load(f)
        except (IOError, json.JSONDecodeError) as e:
            print(f"Error reading JSON file {json_file}: {e}", file=sys.stderr)
            continue
        
        # get genome name (should be first and ONLY key at start of json file)
        if len(list(myk_result.keys())) > 1:
            print(f"Warning: More than one result in mykrobe output file {json_file}, using first result", file=sys.stderr)
        
        genome_name = list(myk_result.keys())[0]
        genome_data = myk_result[genome_name]
        lineage_data = genome_data["phylogenetics"]
        
        lineage_info, spp_call = extract_lineage_info(lineage_data, genome_name, config)
        genome_qrdr_info = extract_qrdr_info(genome_data, genome_name, spp_call, config)
        combined_info = {**lineage_info, **genome_qrdr_info}
        results_data.append(combined_info)

    # Define the column order for output (use species-specific mutations)
    base_columns = [
        "genome", "species", "final genotype", "confidence", 
        "num QRDR"
    ]
    mutation_columns = config['qrdr_mutations']
    coverage_columns = [
        "phylogroup_coverage", "species_coverage", 
        "lowest support for genotype marker", 
        "poorly supported markers", "max support for additional markers", 
        "additional markers", "node support"
    ]
    
    columns = base_columns + mutation_columns + coverage_columns

    # Write results to TSV file
    output_filename = args.prefix + "_predictResults.tsv"
    try:
        with open(output_filename, 'w', newline='') as tsvfile:
            writer = csv.DictWriter(tsvfile, fieldnames=columns, delimiter='\t')
            writer.writeheader()
            for row_data in results_data:
                writer.writerow(row_data)
        print(f"Results written to {output_filename}")
    except IOError as e:
        print(f"Error writing output file {output_filename}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
