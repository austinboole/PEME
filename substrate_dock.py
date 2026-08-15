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


                x = line[31:38].strip() #Grabs the x coordinate of atom
                y = line[39:46].strip() #Grabs the y coordinate of atom
                z = line[47:54].strip() #Grabs the z coordinate of atom

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
                    

REACTION_RULES = {
    "phosphate_hydrolysis":{
        "target_element": "P",
        "required_neighbours": ["O","O","O","O"]
    }


}         

def find_target_atom(atoms,graph,rule):
    target_elem = rule["target_element"]
    req_neighbours = sorted(rule["required_neighbours"])

    matched_targets = []

    for atom_id, data in atoms.items():
        if data["element"] != target_elem:
            continue

        neighbours = graph.get(atom_id,[])
        if len(neighbours) != len(req_neighbours):
            continue

        neighbour_elements = []
        for nbr_id in neighbours:
            if nbr_id in atoms:
                neighbour_elements.append(atoms[nbr_id]["element"])

        if sorted(neighbour_elements) == req_neighbours:
            matched_targets.append(atom_id)

    return matched_targets



# --------------------
#      EXECUTION
# --------------------
file_path = r"C:\Users\Austin Boole\OneDrive\Documents\PEME\inputs\a73c0b85-d39c-4b0d-be5e-26a31ef19321.pdb" #Substrate is G6P in this instance

atoms, graph = parse_pdb(file_path)
targets = find_target_atom(atoms, graph, REACTION_RULES["phosphate_hydrolysis"])
print("Target Atom IDs:",targets)
        
    
