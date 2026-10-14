import os
direc = "/Users/rithwikghanta/Documents/Documents/research/network_polymerization/workspace/"

for job_id in os.listdir(direc):
    jdir = direc + job_id + "/"
    if not os.path.isdir(jdir):
        continue
    if os.path.isdir(jdir + "deform.gsd"):
        continue
    file = jdir + "polymerize_trunc.gsd"
    os.system("python ./scripts/deform.py " + file + " -1")
    break