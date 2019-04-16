import sys
import json
import os

kallistofolder = sys.argv[1]
with open("kallisto_results_filtered.txt","w") as outfile:
    for runID in os.listdir(kallistofolder):
        with open(kallistofolder+"/"+runID + "/run_info.json","r") as run_info_file:
            run_info = json.load(run_info_file)
            p_alignment = run_info["p_pseudoaligned"]
            if p_alignment > 25.0:
                outfile.write(kallistofolder + "/" + runID + "/abundance.h5\n")
