import sys, re, os, gsd, json
import gsd.hoomd
import pandas as pd

#for i in os.listdir()
direc = "/Users/rithwikghanta/Documents/Documents/research/network_polymerization/workspace/"
TE_old = False
failed = []
for job_id in os.listdir(direc):
    # initialize percolation_file dataframe
    if not os.path.isfile(direc + job_id + "/signac_statepoint.json"):
        continue
    #if os.path.isfile(direc + job_id + "/polymerize_trunc.json"):
    #    continue
    
    # the files were truncated for TE_only – take this into account...
    with open(direc + job_id + "/signac_statepoint.json") as f:
            statepoint = json.load(f)
            if (statepoint["chain_side_reaction_probability"]==0):
                TE_old = True
                input_file = direc + job_id + "/polymerize1.gsd"
            else:
                TE_old = False
                input_file = direc + job_id + "/polymerize.gsd"
    print(job_id)

    # read percolation file 
    percolation_file = direc +job_id + "/percolation_data.txt"
    df = pd.read_table(percolation_file, sep=" ", skiprows = None)
    df.dropna(inplace=True, axis=1)
    df.columns = ["frame", "percolation_dimension", "largest_cluster", "molecule_size"]

    # collect the first frame
    frames = [0]
    # find the first instance of percolation for each dimension and record it
    for dim in [1,2,3]:
        first_instance = df[df["percolation_dimension"] == dim]["frame"].min()
        if first_instance > 0:
            frames.append(int(first_instance))
    # collect the last frame
    frames.append(df["frame"].max())
    print(frames)
    # again, fix the frame numbers if these are the truncated polymerize
    if TE_old == True:
        frames = [int(frame//10) for frame in frames]

    # take the corrresponding frames from polymerize and truncate.
    print(frames)
    try:
        output_file = direc + job_id + "/polymerize_trunc.gsd"
        trajectory = gsd.hoomd.open(name=input_file, mode='r')
        print(len(trajectory))
        t =  gsd.hoomd.open(name=output_file, mode='w')
        # Iterate over all frames of the sequence.
        for i in frames:
            t.append(trajectory[int(i)])

        sys.stdout.write("\nDONE\n")
        sys.stdout.flush()
    except:
        print("beep boop :()")
        failed.append(job_id)
        trajectory = gsd.hoomd.open(name=input_file, mode='r')
        t =  gsd.hoomd.open(name=output_file, mode='w')
        for i in [0, len(trajectory)-1]:
            t.append(trajectory[int(i)])
        sys.stdout.write("\nDONE\n")
        sys.stdout.flush()
print(failed)



