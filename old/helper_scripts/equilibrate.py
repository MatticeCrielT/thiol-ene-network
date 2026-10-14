import os
import json

nonequi = []
for sd in os.listdir("./workspace/"):
	if os.path.isfile("./workspace/" + sd + "/equi.gsd"):
		continue
	try:
		os.system("python project.py run -o equilibrate -n 1 -j " + sd)
	except:
		nonequi.append(sd)
		os.system("rm -r " + "./workspace/" + sd + "/equi.gsd")
print(nonequi)	
