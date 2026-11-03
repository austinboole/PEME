import numpy as np

def parse_pdb(file_path):
    atoms = {}
    graph = {}

    with open(file_path,'r') as f:
        lines = f.readlines()

        for line in lines:
            if line.startswith("HETATM"):
                atom_id = int(line[7:11].strip())  #Grabs the serial number of each atom
                atom_name = line[77:78].strip()    #Grabs the element of each atom


                x = float(line[31:38].strip()) #Grabs the x coordinate of atom
                y = float(line[39:46].strip()) #Grabs the y coordinate of atom
                z = float(line[47:54].strip()) #Grabs the z coordinate of atom

                atoms[atom_id]={
                    "element":atom_name,
                    "coord": np.array([x,y,z],dtype=np.float64)
                }

                graph[atom_id] = []

        for line in lines:
            if line.startswith("CONECT"):
                parts = [int(p) for p in line.split()[1:]] #Cuts out the 'CONECT' and converts the list into integers

                source_id = parts[0]
                neighbour_ids = parts[1:]

                for neighbour_id in neighbour_ids:
                    if neighbour_id not in graph[source_id]:      
                        graph[source_id].append(neighbour_id)   #Adds the neighbour atom to the source atom if not already present
                    if source_id not in graph[neighbour_id]:
                        graph[neighbour_id].append(source_id)   #Adds the original source atom to the neighbour atom if not already present

        return atoms,graph
                    

REACTION_RULES = {     #Currently not compatible with methanol, methanonic acid, single carbons etc
    "phosphate_hydrolysis":{
        "target_element": "P",
        "leaving_group": {
            "element": "O",
            "bonded_to": ["C","P"]
        },
        "attack_dist": 2.9
    },
}         

def find_target_atom(atoms,graph,rule):
    target_elem = rule["target_element"]             
    lg_rule = rule["leaving_group"] 

    matched_targets = []

    for atom_id, data in atoms.items():
        if data["element"] != target_elem:      #Skips if element does not match
            continue

        neighbours = graph.get(atom_id,[])

        for nbr_id in neighbours:
            if nbr_id not in atoms:
                continue

            if atoms[nbr_id]["element"] == lg_rule["element"]:

                is_leaving_group = False

                for secondary_id in graph.get(nbr_id,[]):
                    if secondary_id == atom_id or secondary_id not in atoms:
                        continue

                    if atoms[secondary_id]["element"] in lg_rule["bonded_to"]:
                        is_leaving_group = True
                        break

                if is_leaving_group == True:
                    matched_targets.append({
                        "target_id": atom_id,
                        "leaving_group_id": nbr_id
                    })

    return matched_targets


def place_nucleophile_coordinate(atoms, target_info, rule):

    t_id = target_info["target_id"]
    lg_id = target_info["leaving_group_id"]
    attack_dist = rule["attack_dist"]

    target_coord = atoms[t_id]["coord"]
    leaving_group_coord = atoms[lg_id]["coord"]

    direction_vector = target_coord - leaving_group_coord

    vector_norm = np.linalg.norm(direction_vector)
    if vector_norm == 0:
        raise ValueError("Target and leaving group coordinates are identical; cannot determine direction.")
    direction_unit_vector = direction_vector / vector_norm

    nucleophile_coord = target_coord + (direction_unit_vector * attack_dist)

    return nucleophile_coord


def export_nucleophile_pdb(nuc_coord, original_file_path, output_file_path):
    with open(original_file_path, 'r') as f:
        lines = f.readlines()

    newlines = []
    for line in lines:
        if line.startswith("HETATM") or line.startswith("ATOM"):
            newlines.append(line)

    nuc_line = f"HETATM 9999  O   NUC A 999    {float(nuc_coord[0]):8.3f}{float(nuc_coord[1]):8.3f}{float(nuc_coord[2]):8.3f}  1.00 20.00           O\n"
    newlines.append(nuc_line)

    with open(output_file_path, 'w') as f:
        f.writelines(newlines)

    print(f"Nucleophile coordinate exported to {output_file_path}")



# --------------------
#      EXECUTION
# --------------------
file_path = r"C:\Users\Austin Boole\OneDrive\Documents\PEME\inputs\PDBGlucose6Phosphate.pdb" #Substrate is G6P in this instance

atoms, graph = parse_pdb(file_path)
rule = REACTION_RULES["phosphate_hydrolysis"]

targets = find_target_atom(atoms, graph, rule)
print("Target Atom IDs:",targets)
for match in targets:
    nucleophile_coord = place_nucleophile_coordinate(atoms, match, rule)
    print(f"Target Atom ID: {match['target_id']}, Leaving Group Atom ID: {match['leaving_group_id']}, Nucleophile Coordinate: {nucleophile_coord}") 

    output_pdb_path = r"C:\Users\Austin Boole\OneDrive\Documents\PEME\inputs\G6P_with_nucleophile.pdb"
    
    export_nucleophile_pdb(nucleophile_coord, file_path, output_pdb_path)



