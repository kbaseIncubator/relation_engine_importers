#!/usr/bin/env python

# This the pass7 version of create_reaction_similarity_matrix.py to handle new
# versions of reactions.tsv and compounds.tsv (mainly dealing with InChi strings
# being replaced by InchiKeys in the compounts).   Aug 2020
#
#  This reads three files, which are currently harcoded in:
#
#  "compounds.tsv"        - this contains relavent biochemical data for all metabolites 
#                           mentioned in reactions
#  "reactions.tsv"        - 
#  "Unique_ModelSEED_Structures.txt" - new for this update.  InChI strings are no longer
#                        in the compounds.tsv file.  InChIkeys are now in that file.  This
#                        file maps compound ID to InChi string (as well as InChIkey and 
#                        SMARTS string - not currently used, but could be
#

# this writes the matrix to stdout, in the format 
#
#  rxn_id_1 rxn_id_2  sim_score1  sim_score2
#     .        .          .            .
#     .        .          .            .
#     .        .          .            .



import csv
import re

from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import DataStructs


inchi_tab = {}
inchikey_tab  = {}

compounds = {}

#
# load file that maps compound ID to InChI string (and key)
# These are loaded into global inchi_tab and inchikey_tab
#

def load_InChI_strings( filename ):

    global inchi_key    
    global inchikey_tab

    with open( filename ) as uf:
        next( uf )       # skip header 
        for line in uf:
            parts = line.strip().split( '\t')
            id = parts[0]
            if parts[1] == "InChI":
                inchi_tab[id] = parts[5]
                if parts[5][0:6] != 'InChI=':
                    print( "bad inchi [{0}]".format( parts[5] ) )
            elif parts[1] == "InChIKey":
                inchikey_tab[id] = parts[5]

#
# load the metabolite data.  Reactions refer to these compound IDs
# and their InChI strings.  This version converts InChI strings to 
# SMARTS strings for reaction similiarity processing later.  
#
# Note: The file which contains InChI strings also contains SMILES
# strings, which gives the option of using those instead (still converting
# them to SMARTS strings, though)
#
def load_compounds( filename ):

    global inchi_key    
    global inchikey_tab

    comps = {}
    bad_count = 0
    blank_count = 0
    bad_key_count = 0

    with open( filename ) as csv_file:
        csvr = csv.DictReader( csv_file, delimiter='\t' )
        for row in csvr:
            id, inchi_k = ( row['id'], row['inchikey'] )
            if ( id in inchikey_tab.keys() ) and ( id in inchi_tab.keys() ):
                if inchi_k == inchikey_tab.get( id ):
                    try:
                        smarts = Chem.MolToSmarts( Chem.MolFromInchi( inchi_tab.get( id ) ) )
                        comps[id] = smarts
                        # print( "output row {0} {1} {2}".format( id, inchi, smarts ) )
                    except Exception:
                        # print( "input row {0} {1}".format( id, inchi ) )
                        # print( "bizarre", sys.exc_info()[0] )
                        bad_count += 1
                else:
                    print( "bad inchi key for {0} {1} != {2}".format( id, inchi_k, inchikey_tab.get( id ) ) )
                    comps[id] = ""
                    bad_key_count += 1
            else:
                comps[id] = ""
                blank_count += 1

    print( "# bad inputs count: {0}".format( bad_count ) )
    print( "# blank inputs count: {0}".format( blank_count ) )
    print( "# bad key inputs: {0}".format( bad_key_count ) )
    return( comps )


def comp_lookup( comp_id ):
    return( compounds.get(comp_id) )

# This loads the reactions table and creates the reactions table and
# creates the fingerprint strings which are used for the similarity 
# calculations.  The two strings are returned in tables (dicts) rxn and diff_fps

def load_reactions( filename ):

    rxns = {}
    diff_fps = {}
    obsolete_count = 0

    with open(filename) as csv_file:
        csvr = csv.DictReader(csv_file, delimiter='\t')

        # for each reaction

        for row in csvr:
            rxn_id, stoich, is_obsolete = (row['id'], row['stoichiometry'], row['is_obsolete'])
            if int(is_obsolete) > 0:
                obsolete_count = obsolete_count+1
                continue
            # print( "{0}  {1}".format( id, stoich) )
            if stoich:                                  # for now, skip blank stoichiometries (if any)
                left_side_compounds = []
                right_side_compounds = []
                smarts = None

                for cstruct in stoich.split(';'):
                    # print( "   cstruct: {0}".format( cstruct ) )
                    #n, compid, state, x, name = re.findall(r'(?:[^:"]|"(?:\\.|[^"])*")+', cstruct)
                    n, compid, state, x, name = re.findall(r'(?:"(?:\\.|[^"])*"|[^:])+', cstruct)[0:5]
                    # print( "     {0}:    {1} {2} {3} {4}".format( cstruct, n, compid, state, name ) )
                    smarts = comp_lookup(compid)
                    if not smarts or (smarts == ""):
                        smarts = None
                        break
                    copies = int(abs(float(n)))
                    if copies == 0:
                        copies = copies + 1
                    if float(n) < 0:
                        for i in range(0, copies):
                            left_side_compounds.append(smarts)
                    else:
                        for i in range(0, copies):
                            right_side_compounds.append(smarts)

                if smarts is not None:
                    # print( "left" )
                    # pprint( left_side_compounds )
                    # for s in left_side_compounds:
                    #    print( s )
                    # print( "right" )
                    # pprint( right_side_compounds )
                    # for s in right_side_compounds:
                    #    print( s )
                    rxn_string = ".".join(left_side_compounds) + ">>" + \
                                 ".".join(right_side_compounds)
                    # print( "rxn string {0}".format( rxn_string ) )
                    fingerprint = AllChem.CreateStructuralFingerprintForReaction(AllChem.ReactionFromSmarts(rxn_string))
                    # pprint( fingerprint )
                    # pprint( dir( fingerprint ) )
                    # pprint( fingerprint.GetNumBits() )
                    # pprint( fingerprint.ToBinary() )
                    diff_fingerprint = AllChem.CreateDifferenceFingerprintForReaction(
                        AllChem.ReactionFromSmarts(rxn_string))
                    # print( "diff_fingerprint is " )
                    # pprint( diff_fingerprint )
                    # pprint( dir( diff_fingerprint ) )
                    # pprint( diff_fingerprint.GetLength() )
                    # pprint( diff_fingerprint.GetNonzeroElements() )
                    # b = diff_fingerprint.ToBinary()
                    # print( type(b) )
                    # pprint( b )
                    rxns[rxn_id] = fingerprint
                    diff_fps[rxn_id] = diff_fingerprint

    print("# obsolete_count = {0}".format(obsolete_count))

    return(rxns, diff_fps)


                              ################
                              # Main Program #
                              ################

# load InChi strings file

load_InChI_strings( "Unique_ModelSEED_Structures.txt" )

#for i in inchi_tab.keys():
#    print( "  {0}   {1}    {2}".format( i, inchikey_tab[i], inchi_tab[i] ) )


# next load compounds and convert to SMARTS and put in table


compounds = load_compounds("compounds.tsv")

# pprint( compounds )

# Next, load reactions, capture reaction strings and replace compound ids with SMARTS

reactions, diffs = load_reactions("reactions.tsv")

rxn_list = list(reactions.keys())   # list() required for python 3
num_rxns = len(rxn_list)
# num_rxns = 10000

#
# Main loop - every pair (these are symmetric scores presumably, so 
# n^2 / 2 pairs need to be considered
#

for i in range(0, num_rxns-1):
    for j in range(i+1, num_rxns):
        rxn_a = rxn_list[i]
        rxn_b = rxn_list[j]
        print("{0} {1} {2} {3}".format(rxn_a, rxn_b,
                                       DataStructs.FingerprintSimilarity(reactions[rxn_a],
                                                                         reactions[rxn_b]),
                                       DataStructs.cDataStructs.TanimotoSimilarity(diffs[rxn_a],
                                                                                   diffs[rxn_b])
                                       )
             )
