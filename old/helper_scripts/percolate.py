import os 
import signac 

# our signac project
project = signac.get_project()
l = ["56fb642c37467c91ec211c993a54b0e0",
"1c6da38a4491c2d8cd4c2768180d0aef",
"665f7582c0248e3657ea1b0a25e47be2",
"0214851522083aa91b4c02c4d2a55495"]

for job in l:
    polymerize_file = "./workspace/" + job+ "/polymerize.gsd"
    os.system("./percolation-analyzer/bin/sample_graph -i " + polymerize_file + " -o " + "./" + job + "_percolation_data.txt")