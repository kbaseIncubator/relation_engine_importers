#!/usr/bin/env python3

# reads compounds.tsv and creates json output for uploading

import sys
import csv
import json

if len( sys.argv ) < 2:
    print( "please provide compound tsv file name" )
    quit()

csvfile = sys.argv[1]


with open( csvfile ) as cf:
    csv = csv.DictReader( cf, delimiter='\t' )
    for rowd in csv:
        print( json.dumps( rowd ) )
        


