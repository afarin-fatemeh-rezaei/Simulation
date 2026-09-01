import json
import glob
import re

json_files = glob.glob("bicrystal_scan_25deg_best.json")

for json_file in json_files:
    match = re.search(r'(\d+)deg', json_file)
    if not match:
        continue
    theta = match.group(1)
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    dat_file = f"atoms_{theta}deg.dat"
    with open(dat_file, 'w') as f:
        for atom in data:
            f.write(f"{atom['x']:.10f} {atom['y']:.10f}\n")
    
    print(f"Converted {json_file} -> {dat_file}")