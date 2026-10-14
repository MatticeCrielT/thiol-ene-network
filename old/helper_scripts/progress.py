import os
import json

reacting_job_ids_dict = {}
reacted_job_ids_dict = {}
no_reaction = []
for job_id in os.listdir("./workspace/"):
    job_path = "./workspace/" + job_id + "/"
    if os.path.isfile(job_path + "polymerize.gsd"):
        with open(job_path + "signac_job_document.json") as job_doc:
            data = json.load(job_doc)
            rm = data["reacted_monomers"]
        if rm < 0.9:
            reacting_job_ids_dict[job_id] = rm
        else:
            reacted_job_ids_dict[job_id] = rm
    elif os.path.isfile(job_path + "polymerize1.gsd"):
        with open(job_path + "signac_job_document.json") as job_doc:
            data = json.load(job_doc)
            rm = data["reacted_monomers"]
        if rm < 0.9:
            reacting_job_ids_dict[job_id] = rm
        else:
            reacted_job_ids_dict[job_id] = rm
    else: 
        no_reaction.append(job_id)
with open("./job_progress/reacting_job_ids.json", 'w') as f:
    json.dump(reacting_job_ids_dict, f)
with open("./job_progress/reacted_job_ids.json", 'w') as f:
    json.dump(reacted_job_ids_dict, f)

print("reacted")
print(reacted_job_ids_dict)

print("reacting")
print(reacting_job_ids_dict)

print("nothing")
print(no_reaction)
