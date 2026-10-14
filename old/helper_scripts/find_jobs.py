import os, json

for job_id in os.listdir("./workspace/"):
    job_path = "./workspace/" + job_id + "/"
    if os.path.isfile(job_path + "polymerize.gsd"):
        with open(job_path + "signac_statepoint.json") as f:
            data = json.load(f)
            if data["chain_side_reaction_probability"] == 0 and data["crosslinker_percent"]==12.5:
                print(job_id)
            else:
                continue
        with open(job_path + "signac_job_document.json") as job_doc:
            data = json.load(job_doc)
            rm = data["reacted_monomers"]
            print(rm)