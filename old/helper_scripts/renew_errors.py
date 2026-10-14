import os

i = 0
while i < 500:
    os.system("python project.py run -o polymerize  -j 29c0424f03e0d3f2257a9049244f7562 -t 90")
    os.remove("./workspace/29c0424f03e0d3f2257a9049244f7562/polymerize.gsd")
    os.system("python project.py run -o polymerize  -j cf4868bca033dcd67548b1f83ed6fa92 -t 90")
    os.remove("./workspace/cf4868bca033dcd67548b1f83ed6fa92/polymerize.gsd")
    i +=1