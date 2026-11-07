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
                    

REACTION_RULES = {
    "phosphate_hydrolysis": {               #Will break in the case of Diesters, Nucleic Acids, and bridging Phosphates, but will work for Monoesters
        "target_element": "P",
        "leaving_group": {
            "element": "O",
            "bonded_to": ["C", "P"]
        },
        "geometric_op": "apply_trigonal_bipyramidal_shift", 
        "pauling_params": {
            "d0": 1.60,      
            "c": 0.60,       
            "bond_order": 0.5 
        }
    }
}
         

def find_target_atom(atoms,graph,rule):
    target_elem = rule["target_element"]             
    lg_rule = rule["leaving_group"] 

    matched_targets = []

    for atom_id, data in atoms.items():
        if data["element"] != target_elem:      #Skips if element does not match
            continue

        neighbours = graph.get(atom_id,[])
        leaving_group_id = None
        coordinating_neighbours = []

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

                if is_leaving_group:
                    leaving_group_id = nbr_id
                    continue

            coordinating_neighbours.append(nbr_id)

        if leaving_group_id is not None:
            matched_targets.append({
                "target_atom_id": atom_id,
                "leaving_group_id": leaving_group_id,
                "coordinating_neighbour_ids": coordinating_neighbours
            })

    return matched_targets


def calc_linus_pauling_tranistion_state_distance(params):
    """
    Calculates the stretched bond length at the transition state using 
    Linus Pauling's bond-order/bond-length relationship: d(n) = d0 - c * ln(n)
    """
    d0 = params["d0"]
    c = params["c"]
    bond_order = params["bond_order"]

    d_n = d0 - c * np.log(bond_order)
    return d_n

def apply_trigonal_bipyramidal_shift(atoms, target_id, leaving_group_id, coordinating_neighbour_ids, d_ts):
    target_coord = atoms[target_id]["coord"]
    lg_coord = atoms[leaving_group_id]["coord"]

    # Establish the vector from the target atom to the leaving group

    apical_vector = target_coord - lg_coord
    vector_norm = np.linalg.norm(apical_vector)  # Normalize
    apical_unit = apical_vector / vector_norm

    updated_coords = {}

    updated_coords[leaving_group_id] = lg_coord
    
    new_target_coord = lg_coord + (apical_unit*d_ts)
    updated_coords[target_id] = new_target_coord

    updated_coords["nucleophile"] = new_target_coord + (apical_unit*d_ts)

    # Calculate the new positions for the coordinating neighbours in the equatorial plane

    for nbr_id in coordinating_neighbour_ids:
        nbr_coord = atoms[nbr_id]["coord"]

        r_vector = nbr_coord - target_coord
        orig_bond_length = np.linalg.norm(r_vector)
        normal_component = np.dot(r_vector, apical_unit) * apical_unit
        in_plane_vector = r_vector - normal_component

        in_plane_norm = np.linalg.norm(in_plane_vector)
        if in_plane_norm > 0:
            in_plane_unit = in_plane_vector / in_plane_norm
            new_nbr_coord = new_target_coord + (in_plane_unit*orig_bond_length)
            updated_coords[nbr_id] = new_nbr_coord
            
            # --- RODRIGUES' ROTATION FOR ATTACHED HYDROGENS ---

            old_bond_vector = nbr_coord - target_coord
            new_bond_vector = new_nbr_coord - new_target_coord

            old_bond_unit = old_bond_vector / np.linalg.norm(old_bond_vector)
            new_bond_unit = new_bond_vector / np.linalg.norm(new_bond_vector)

            rotation_axis = np.cross(old_bond_unit, new_bond_unit)
            axis_norm = np.linalg.norm(rotation_axis)

            for attached_id in graph.get(nbr_id, []):
                if attached_id in atoms and atoms[attached_id]["element"] == "H":
                    attached_vector = atoms[attached_id]["coord"] - nbr_coord
                    if axis_norm > 1e-6:
                        k = rotation_axis / axis_norm
                        cos_theta = np.dot(old_bond_unit, new_bond_unit)
                        theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))
                        sin_theta = np.sin(theta)

                        rotated_vector = (attached_vector * cos_theta + np.cross(k, attached_vector) * sin_theta + k * np.dot(k, attached_vector) * (1 - cos_theta))
                        updated_coords[attached_id] = new_nbr_coord + rotated_vector
                    else:
                        updated_coords[attached_id] = new_nbr_coord + attached_vector

        else:
            raise ValueError(f"Atom {nbr_id} is collinear with the apical axis.")
    return updated_coords


def export_transition_state_pdb(orig_file_path, new_file_path, updated_coords):
    with open(orig_file_path, 'r') as f:
        lines = f.readlines()

        output_lines = []
        nuc_added = False  # Flag to track if nucleophile has been added
        for line in lines:

            if line.startswith("CONECT") and "nucleophile" in updated_coords and nuc_added == False:
                nuc_coord = updated_coords["nucleophile"]
                output_lines.append(f"HETATM 9999  O   NUC A 999    {nuc_coord[0]:8.3f}{nuc_coord[1]:8.3f}{nuc_coord[2]:8.3f}  1.00 20.00           O\n")
                nuc_added = True  # Set the flag to True after adding the nucleophile

            if line.startswith("HETATM") or line.startswith("ATOM"):
                atom_id = int(line[7:11].strip())
                if atom_id in updated_coords:
                    coord = updated_coords[atom_id]
                    line = f"{line[:30]}{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}{line[54:]}"

            output_lines.append(line)

    with open(new_file_path, 'w') as f:
        f.writelines(output_lines)

    print(f"Transition state PDB exported to {new_file_path}")




# --------------------
#      EXECUTION
# --------------------

file_path = r"C:\Users\Austin Boole\OneDrive\Documents\PEME\inputs\PDBGlucose6Phosphate.pdb" #Substrate is G6P in this instance
output_path = r"C:\Users\Austin Boole\OneDrive\Documents\PEME\inputs\TransitionState_G6P.pdb"

atoms, graph = parse_pdb(file_path)
rule = REACTION_RULES["phosphate_hydrolysis"]

targets = find_target_atom(atoms, graph, rule)
if targets:
    print("Target Atom IDs:",targets)

    t = targets[0]
    d_ts = calc_linus_pauling_tranistion_state_distance(rule["pauling_params"])

    updated_coords = apply_trigonal_bipyramidal_shift(
        atoms = atoms,
        target_id = t["target_atom_id"],
        leaving_group_id = t["leaving_group_id"],
        coordinating_neighbour_ids = t["coordinating_neighbour_ids"],
        d_ts = d_ts
    )

    print("Updated Coordinates:")
    for atom_id, coord in updated_coords.items():
        if atom_id == "nucleophile":
            print(f"Nucleophile Position: {coord}")
        else:
            print(f"Atom ID {atom_id} ({atoms[atom_id]['element']}): {coord}")


    export_transition_state_pdb(file_path, output_path, updated_coords)



else:
    print("No suitable target atoms found for the specified reaction rule.")




