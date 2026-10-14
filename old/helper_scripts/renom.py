import os 

direc = "./workspace/"
for job_id in os.listdir(direc):
    jdir = direc + job_id + "/"
    if not os.path.isdir(jdir):
        continue
    for f in os.listdir(jdir):
        if "polymerize1" in f:
            os.system("mv " + jdir + f + " " + jdir + "polymerize.gsd")
