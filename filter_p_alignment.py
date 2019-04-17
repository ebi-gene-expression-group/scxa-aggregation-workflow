import sys
import json
import os

# Filter out the cells for which less than minimum_p % of the reads could be pseudo aligned to the transcriptome by kallisto.

if len(sys.argv)<3:
    print("Not enough arguments. Usage: python filter_p_alignment.py kallistofolder minimum_p")
kallistofolder = sys.argv[1]
minimum_p = sys.argv[2]
with open("kallisto_results_filtered.txt","w") as outfile:
    for runID in os.listdir(kallistofolder):
        with open(kallistofolder+"/"+runID + "/run_info.json","r") as run_info_file:
            run_info = json.load(run_info_file)
            p_alignment = run_info["p_pseudoaligned"]
            if p_alignment > minimum_p:
                outfile.write(kallistofolder + "/" + runID + "/abundance.h5\n")
