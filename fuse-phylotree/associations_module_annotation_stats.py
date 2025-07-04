#!/bin/python3

import argparse
from pathlib import Path
import csv
from ete3 import Tree
from scipy.stats import fisher_exact


#==============================================================================
# Read and load files
#==============================================================================

import csv

def _split_modules(cell: str) -> list[str]:
    """
    Turn 'B10:1.0|B12:0.9' → ['B10', 'B12'].
    Works for empty strings too.
    """
    cell = cell.strip()
    if not cell:
        return []
    return [part.split(":", 1)[0] for part in cell.split("|") if part]

def load_evo_mod_func(filepath: str) -> dict:
    """
    Parse the evolution table with columns:
        gene, modules_present, function_present,
        modules_gained,  function_gained,
        modules_lost,    function_lost
    Returns {gene: {...}} where every *modules_* entry is a list of names
    (frequencies discarded).
    """
    dict_gene_evoModFunc = {}

    with open(filepath, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)

        for row in reader:
            gene = row["gene"].strip()

            dict_gene_evoModFunc[gene] = {
                "modules_present":  _split_modules(row["modules_present"]),
                "function_present": row["function_present"].strip().split("|")
                                    if row["function_present"].strip() else [],

                "module_gained":    _split_modules(row["modules_gained"]),
                "function_gained":  row["function_gained"].strip().split("|")
                                    if row["function_gained"].strip() else [],

                "module_lost":      _split_modules(row["modules_lost"]),
                "function_lost":    row["function_lost"].strip().split("|")
                                    if row["function_lost"].strip() else [],
            }

    return dict_gene_evoModFunc


#==============================================================================
# Compute stats
#==============================================================================

def compute_module_function_stats(dict_gene_evoModFunc: dict) -> list:

    # Collect all unique modules and functions
    all_modules = set()
    all_functions = set()
    for gene_data in dict_gene_evoModFunc.values():
        all_modules.update(gene_data['modules_present'])
        all_functions.update(gene_data['function_present'])

    results = []

    # Prepare presence mappings for modules and functions
    module_to_genes = {mod: set() for mod in all_modules}
    function_to_genes = {func: set() for func in all_functions}
    module_gain_map = {mod: set() for mod in all_modules}
    function_gain_map = {func: set() for func in all_functions}

    for gene, data in dict_gene_evoModFunc.items():
        for mod in data['modules_present']:
            module_to_genes[mod].add(gene)
        for func in data['function_present']:
            function_to_genes[func].add(gene)
        for mod in data['module_gained']:
            module_gain_map[mod].add(gene)
        for func in data['function_gained']:
            function_gain_map[func].add(gene)

    all_genes = set(dict_gene_evoModFunc.keys())

    # Compute stats for each module-function pair
    for mod in all_modules:
        for func in all_functions:
            genes_with_mod = module_to_genes[mod]
            genes_with_func = function_to_genes[func]

            # Jaccard index
            intersection = genes_with_mod & genes_with_func
            union = genes_with_mod | genes_with_func
            jaccard_index = len(intersection) / len(union) if union else 0.0

            # Fisher's Exact Test
            a = len(intersection)  # both mod and func
            b = len(genes_with_mod - genes_with_func)  # mod only
            c = len(genes_with_func - genes_with_mod)  # func only
            d = len(all_genes - (genes_with_mod | genes_with_func))  # neither
            contingency_table = [[a, b], [c, d]]
            try:
                oddsratio, p_value = fisher_exact(contingency_table)
            except:
                oddsratio, p_value = float('nan'), 1.0

            # Co-emergence
            co_emergence = bool(module_gain_map[mod] & function_gain_map[func])

            # Only keep non-null associations 
            if jaccard_index > 0.00:
                results.append([mod, func, co_emergence, jaccard_index, p_value, oddsratio])

    return results


#==============================================================================
# Write output
#==============================================================================

def write_module_function_stats_to_csv(results: list, output_file: str) -> None:
    headers = ['module', 'function', 'co_emergence', 'co-presence_jaccard_index', 'co-presence_FET_p_value', 'co-presence_FET-odds_ratio']  
    with open(output_file, mode='w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        for row in results:
            writer.writerow(row)

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--modules_functions_evolution",
                        help = "The final csv file containing the modules/functions presences and gained for the different gene nodes (eg., 1_modules_and_functions_evolution.csv)",
                        type=str)
    # Not used for now - simply looking at co-presences
    #parser.add_argument("--gene_tree",
    #                    help = "The final gene tree corresponding the csv files with modules/functions presences and gained (eg., 0_gene.tree)",
    #                    type=str) 
    parser.add_argument("--output",
                        help = "Custom output csv filename (only associations with Jaccard index > 0.00)",
                        type=str,
                        default='X_module_function_assoc_stats.csv')
    # Filter ?
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    # Get arguments
    args = parser()
    fn_evo_mod_func = Path(args.modules_functions_evolution).resolve()
    fn_output = Path(args.output).resolve()
    #fn_gene_tree = Path(args.gene_tree).resolve()
    # Load data
    dict_gene_evoModFunc = load_evo_mod_func(fn_evo_mod_func)
    #gene_tree = Tree(fn_gene_tree, format=1)
    # Compute associations stats
    results = compute_module_function_stats(dict_gene_evoModFunc)
    # Write output
    write_module_function_stats_to_csv(results, str(fn_output))

if __name__ == '__main__':
    main()
