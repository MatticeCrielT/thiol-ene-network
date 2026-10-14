import os 

for i in os.listdir("./workspace/"):
    if os.path.isdir("./workspace/" + i) and os.path.isfile("./workspace/" + i + "/alternating_q.json"):
        print("hi")
        os.rename("./workspace/" + i + "/equi.gsd", "./workspace/" + i + "/equi_1.gsd")