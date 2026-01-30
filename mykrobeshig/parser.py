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
import re
import sys
import csv
from importlib import resources
from argparse import ArgumentParser


def get_arguments():
    parser = ArgumentParser(description='Parse mykrobe predict JSON files for Shigella sonnei or Shigella flexneri')

    parser.add_argument('--jsons', required=True, nargs='+', help='JSON files output from mykrobe predict')
    parser.add_argument('--prefix', required=True, help='prefix for output files')
    parser.add_argument('--force', required=False, action='store_true', help='If you have used --force with Mykrobe to get a lineage call despite species mapping failing, use that same logic here. Lineage will be reported even though species failed.')
    return parser.parse_args()


def get_paramters_for_spp(spp_call, scheme):
    """
    Determine the appropriate parameters based on the species call from Mykrobe.
    If it's not one of these species, return None (will be processed as "Unknown")
    """
    # Map Mykrobe species names to config keys
    species_config = {
    'Shigella_flexneri': {
        'scheme': 'flexneri',
        'display_name': 'S. flexneri',
        'coverage_threshold': 85.5,
        'phylo_group': 'Ecoli_Shigella',
        'qrdr_mutations': ["parC_S80I", "gyrA_S83L", "gyrA_S83A", "gyrA_D87G", "gyrA_D87N", "gyrA_D87Y"]
    },
    'Shigella_sonnei': {
        'scheme': 'sonnei',
        'display_name': 'S. sonnei',
        'coverage_threshold': 90,
        'phylo_group': 'Ecoli_Shigella',
        'qrdr_mutations': ["parC_S80I", "gyrA_S83L", "gyrA_S83A", "gyrA_D87G", "gyrA_D87N", "gyrA_D87Y"]
    }
}
    if spp_call in species_config.keys():
        return species_config[spp_call]
    else:
        # If unknown or unrecognized species
        # return the value for the scheme
        for key, value in species_config.items():
            if value['scheme'] == scheme:
                return value
        #return None


def determine_scheme(genome_data):

    """Determine which genotyping scheme was used based on the probe sets. Returns the species name of the scheme, or None if undetermined."""
    probe_sets = genome_data["probe_sets"]
    flex_pattern = re.compile(r'flexneri\..+\.fa\.gz$')
    sonnei_pattern = re.compile(r'sonnei\..+\.fa\.gz$')

    flex_count = sum(1 for f in probe_sets if flex_pattern.search(f))
    sonnei_count = sum(1 for f in probe_sets if sonnei_pattern.search(f))
    
    if sonnei_count >= 1:
        return "sonnei"
    elif flex_count >= 1:
        return "flexneri"
    else:
        return None
    

def extract_qrdr_info(genome_data, genome_name, spp_call, config):
    """Extract QRDR mutation information for the genome"""
    qrdr_out_dict = {}
    if config:
        qrdr_possible = config['qrdr_mutations']
    else:
        qrdr_possible = ["parC_S80I", "gyrA_S83L", "gyrA_S83A", "gyrA_D87G", "gyrA_D87N", "gyrA_D87Y"]
        qrdr_out_dict['genome'] = genome_name
        for allele in qrdr_possible:
            qrdr_out_dict[allele] = 'NA'
        qrdr_out_dict['num QRDR'] = 'NA'
        return qrdr_out_dict

    # can only extract qrdr info if sample matches the species, otherwise set 'NA' for all calls
    #if spp_call == config['species_name']:
    try:
        qrdr_data = genome_data["susceptibility"]["ciprofloxacin"]
    except KeyError:
        qrdr_out_dict['genome'] = genome_name
        for allele in qrdr_possible:
            qrdr_out_dict[allele] = 'NA'
        qrdr_out_dict['num QRDR'] = 'NA'
        return qrdr_out_dict
    #else:
    #    qrdr_out_dict['genome'] = genome_name
    #    for allele in qrdr_possible:
    #        qrdr_out_dict[allele] = 'NA'
    #    qrdr_out_dict['num QRDR'] = 'NA'
    #    return qrdr_out_dict

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


def inspect_calls(full_lineage_data, spp):
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
        # regardless of the call, get info 
        call_details = full_lineage_data['calls'][best_genotype][level]
        # check that there is something there
        if call_details:
            # need to do this weird thing to grab the info without knowing the key name
            call_details = call_details[list(call_details.keys())[0]]
            ref = call_details['info']['coverage']['reference']['median_depth']
            alt = call_details['info']['coverage']['alternate']['median_depth']
            # calculate percent support for marker
            try:
                percent_support = alt / (alt + ref)
            except ZeroDivisionError:
                percent_support = 0
            
            marker_string = level + ' (' + str(best_calls[level]) + '; ' + str(alt) + '/' + str(ref) + ')'
            final_markers.append(marker_string)
        # if the value is null, just report 0 (indicates that no SNV detected, either ref or alt?)
        else:
            # note we do not have a markers for lineage5.1 in sonnei so don't report this as 0
            if level != 'lineage5.1' and spp == 'Shigella_sonnei':
                lowest_within_genotype_percents[0] = level
                marker_string = level + ' (0)'
                poorly_supported_markers.append(marker_string)
                final_markers.append(marker_string)
        # if call is 1 then that is fine, count to determine confidence later
        # if call is 0.5, then get info
        # if call is 0, there will be no info in the calls section, so just report 0s everywhere
        if best_calls[level] < 1 or (best_calls[level] < 1 and level != 'lineage5.1' and spp == 'Shigella_sonnei'):
            # then it must be a 0 or a 0.5
            # report the value (0/0.5), and also the depth compared to the reference
            lowest_within_genotype_percents[percent_support] = level
            poorly_supported_markers.append(marker_string)

    # determining final confidence is based ONLY on the actual genotype, not incongruent genotype calls
    # if there is a lineage5.1 in the call, then we need to remove one of the 0s in best_calls_vals
    if 'lineage5.1' in best_calls.keys() and spp == 'Shigella_sonnei':
        best_calls_vals.remove(0)
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

    # make a list of all possible quality issues (incongruent markers, or not confident calls within the best geno)
    non_matching_markers = []
    non_matching_supports = []
    # we now want to report any additional markers that aren't congruent with our best genotype
    #ie if 3.6.1 is the best genotype, but we also have a 3.7.29 call, we need to report the 3.7 and 3.29 markers as incongruent
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


def phylogroup_no_match(phylo_grp_call, phylo_grp_data, genome_name, scheme):
    """Create output dictionary for genomes that do not match the Ecoli/Shigella phylogroup"""
    if phylo_grp_call != 'Unknown':
        phylo_grp_percentage = phylo_grp_data[phylo_grp_call]['percent_coverage']
        species_val = 'Phylogroup ' + phylo_grp_call
    else:
        phylo_grp_percentage = 'NA'
        species_val = 'Not E. coli/Shigella'
    
    # create the output dictionary
    out_dict = {
        'genome': genome_name, 
        'species': species_val,
        'final genotype': 'NA',
        'name': 'NA',
        'scheme': scheme,
        'confidence': 'NA',
        'phylogroup_coverage': phylo_grp_percentage,
        'species_coverage': 'NA',
        'lowest support for genotype marker': '', 
        'poorly supported markers': '', 
        'node support': '', 
        'max support for additional markers': '', 
        'additional markers': ''
    }
    return out_dict


def extract_lineage_info(lineage_data, genome_name, config, sonnei_name_dict, force=False):
    """Extract lineage and species information"""
    
    spp_data = lineage_data['species']
    spp_call = list(spp_data.keys())[0]
    phylo_grp_data = lineage_data['phylo_group']
    
    # grab the percentages for species and phylogroup
    if not force:
        spp_percentage = spp_data[spp_call]["percent_coverage"]
        phylo_grp_percentage = phylo_grp_data['Ecoli_Shigella']['percent_coverage']
    # if we're forcing, set these to -1
    else:
        spp_percentage = -1
        phylo_grp_percentage = -1

    
    # if the percentage is below threshold, then don't parse any further (if we aren't forcing a call)
    # but we call this as 'below threshold for species', as it is
    # Ecoli/Shigella, just not a high enough match to the species we expect
    if spp_percentage < config['coverage_threshold'] and not force:
        out_dict = {
            'genome': genome_name, 
            'species': 'Unknown E. coli/Shigella', 
            'final genotype': 'NA',
            'name': 'NA',
            'scheme': config['scheme'],
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

    # Otherwise continue if it's one of sonnei or flexneri
    # also continue if we're using 'force', even if below threshold
    lineage_out_dict = {'genome': genome_name}

    try:
        genotype_calls = lineage_data['lineage']['lineage']
    except KeyError:
        genotype_calls = []
    
    # if there are no calls, populate with none
    if len(genotype_calls) == 0:
        lineage_out_dict['final genotype'] = 'uncalled'
        lineage_out_dict['name'] = 'NA'
        lineage_out_dict['confidence'] = 'NA'
        lineage_out_dict['lowest support for genotype marker'] = ''
        lineage_out_dict['poorly supported markers'] = ''
        lineage_out_dict['max support for additional markers'] = ''
        lineage_out_dict['additional markers'] = ''
        lineage_out_dict['node support'] = ''
    else:
        # inspect calls for genotype(s)
        best_genotype, confidence, lowest_support_val, poorly_supported_markers, non_matching_support, non_matching_markers, final_markers = inspect_calls(lineage_data['lineage'], spp_call)
        lineage_out_dict['final genotype'] = best_genotype
        lineage_out_dict['confidence'] = confidence
        lineage_out_dict['lowest support for genotype marker'] = lowest_support_val
        lineage_out_dict['poorly supported markers'] = '; '.join(poorly_supported_markers)
        lineage_out_dict['max support for additional markers'] = non_matching_support
        lineage_out_dict['additional markers'] = '; '.join(non_matching_markers)
        lineage_out_dict['node support'] = '; '.join(final_markers)
        if spp_call == "Shigella_sonnei":
            try:
                lineage_out_dict['name'] = sonnei_name_dict[best_genotype]
            except KeyError:
                lineage_out_dict['name'] = '-'
        else:
            lineage_out_dict['name'] = '-'
    
    # add conserved info
    # if forcing, indicate this
    if force:
        lineage_out_dict['species'] = 'forced call'
    else:
        lineage_out_dict['species'] = config['display_name']
    lineage_out_dict['scheme'] = config['scheme']
    lineage_out_dict['phylogroup_coverage'] = phylo_grp_percentage
    lineage_out_dict['species_coverage'] = spp_percentage

    return lineage_out_dict, spp_call


def main():
    args = get_arguments()

    # create sonnei_name_dict, key=mykrobe lineage name, value=sonnei human readable name
    # only used for sonnei genotypes
    sonnei_name_dict = {}
    with resources.open_text('mykrobeshig.data', 'alleles_sonnei.txt') as lineage_names:
        for line in lineage_names:
            fields = line.strip().split('\t')
            sonnei_name_dict[fields[4]] = fields[3]
    
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
        
        # Extract what species the genome has been called against
        # by looking at the probe_sets values - if 'sonnei' then has been called
        # against sonnei scheme, if 'flexneri' then flex scheme, if neither then
        # this parser script shouldn't apply
        scheme = determine_scheme(genome_data)
        if not scheme:
            # if we can't determine the scheme, then skip to next genome
            print(f"Warning: Genotyping scheme for genome {genome_name} is not flexneri or sonnei, skipping.", file=sys.stderr)
            continue
        
        # Now that we know we're typing either a supposed sonnei or flexneri, extract the phylogroup and species calls Mykrobe made
        spp_data = lineage_data['species']
        spp_call = list(spp_data.keys())[0]
        phylo_grp_data = lineage_data['phylo_group']
        phylo_grp_call = list(phylo_grp_data.keys())[0]

        # if our phylogroup isn't Ecoli/Shigella, then skip extracting any lineage info
        # only do this if we're not 'forcing' a call
        if phylo_grp_call != 'Ecoli_Shigella' and not args.force:
            out_dict = phylogroup_no_match(phylo_grp_call, phylo_grp_data, genome_name, scheme)
            qrdr_no_call = extract_qrdr_info(genome_data, genome_name, spp_call, None)
            out_dict.update(qrdr_no_call)
            results_data.append(out_dict)
            continue
        else:
            # Grab the correct parameters for this sceheme
            spp_config = get_paramters_for_spp(spp_call, scheme)
            if spp_config is None and args.force:
                # if we're forcing, then set config to whatever the scheme is
                if scheme == 'sonnei':
                    spp_config = get_paramters_for_spp('Shigella_sonnei')
                elif scheme == 'flexneri':
                    spp_config = get_paramters_for_spp('Shigella_flexneri')
                else:
                    raise
            
            lineage_info, spp_call_returned = extract_lineage_info(lineage_data, genome_name, spp_config, sonnei_name_dict, args.force)
            genome_qrdr_info = extract_qrdr_info(genome_data, genome_name, spp_call_returned, spp_config)
            combined_info = {**lineage_info, **genome_qrdr_info}
            results_data.append(combined_info)

    # Define the columns
    base_columns = [
        "genome", "species", "final genotype", "name", "scheme", "confidence", 
        "num QRDR"
    ]

    qrdr_cols = ["parC_S80I", "gyrA_S83L", "gyrA_S83A", "gyrA_D87G", "gyrA_D87N", "gyrA_D87Y"]

    coverage_columns = [
        "phylogroup_coverage", "species_coverage", 
        "lowest support for genotype marker", 
        "poorly supported markers", "max support for additional markers", 
        "additional markers", "node support"
    ]
    
    # define the column order for the output
    columns = base_columns + qrdr_cols + coverage_columns

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
