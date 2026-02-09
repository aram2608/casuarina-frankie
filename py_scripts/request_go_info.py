import requests
from requests import Response
import json
import time

endpoint = "https://api.geneontology.org/api/go"

with open("go_counts.json", "r") as f:
    obj = f.read()

obj = json.loads(obj)

out = {}

for k in obj.keys():
    res: Response = requests.get(f"{endpoint}/{k}")
    out[k] = res.json()
    time.sleep(0.5)

out = json.dumps(out, indent=4)

with open("data.json", "w") as f:
    f.write(out)