import sys
from os import path

r1_path = sys.argv[3]
# R2 is optional: when absent, the manifest is written for a single-end sample
# (forward read only).
r2_path = sys.argv[4] if len(sys.argv) > 4 else None

#parts = path.basename(r1_path).replace('.gz', '').replace('.fq', '').replace('.fastq', '').split('_')
#to_use_in_name = [x for x in parts if not x.lstrip('SLR').isnumeric()]
sample_name = path.basename(r1_path).split('_')[0].replace('.gz', '').replace('.fq', '').replace('.fastq', '')

manH = open(sys.argv[1], 'w')
manH.write( "sample-id,absolute-filepath,direction\n" )
manH.write( f"{ sample_name },{ r1_path },forward\n" )
if r2_path:
    manH.write( f"{ sample_name },{ r2_path },reverse\n" )
manH.close()

metaH = open(sys.argv[2], 'w')
metaH.write( "sample-id\tcondition\n" )
metaH.write( f"{ sample_name }\t{ sample_name }\n" )
metaH.close()