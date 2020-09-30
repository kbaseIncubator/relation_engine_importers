"""
Generate json import file for the master compound list

Can be run with: python generate_compound.py compound.tsv path/to/output/dir

Get compounds.tsv from: https://github.com/kbaseapps/BiochemistryAPI/blob/master/data/compounds.tsv
"""
import os
import sys
import csv

from utils.write_import_file import write_import_file

_compound_vert_name = "rxn_compound"

# Transformations on the header names
# Eg. given a header "id", save it as a field named "_key" in the db
header_transforms = {
    'id': '_key',
    'kegg id': 'kegg_id',
    'ms id': 'ms_id'
}





def iterate_tsv_rows(file_path):
    """
    Iterate over every row in a TSV file and import each as a compound document.
    For each compound, we create:
    - A compound
    - A set of edges from the compound to all constituent compounds
    - A set of gene complexes based on the 'gpr' boolean expression of required genes
    - Edges from the gene complexes to the genes
    - Edges from the gene complexes to the compound that gets produced
    """
    with open(file_path, newline="") as csv_fd:
        reader = csv.reader(csv_fd, delimiter="\t", quotechar='"')
        headers = [header_transforms.get(h, h) for h in next(reader)]
        for row in reader:
            compound_gen = gen_compound(row, headers)
            for result in compound_gen:
                yield result


def gen_compound(row, headers):
    """Import a single compound from a row in a TSV file."""
    compound = {}  # type: dict
    for (idx, col) in enumerate(row):
        compound[headers[idx]] = col
    if "_key" not in compound:
        return
    compound['id'] = compound['_key']
    yield compound


def _fatal(msg):
    """Write error message to standard out and exit 1."""
    sys.stderr.write(msg + "\n")
    sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Pass in two arguments:")
        print("  - path to compounds.tsv file")
        print("  - path to output directory")
        _fatal("Invalid arguments")
    tsv_path = sys.argv[1]
    output_dir = sys.argv[2]
    if not os.path.exists(tsv_path):
        _fatal("Path does not exist: %s" % tsv_path)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, _compound_vert_name + ".json")
    write_import_file(iterate_tsv_rows(tsv_path), output_file)
    print("done.")
