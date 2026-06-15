import numpy as np
import matplotlib.pyplot as plt

class ProteinBackboneGenerator:
    def __init__(self,radius=2.3,rise_p_residue=1.5,residues_p_turn=3.6):
        self.r = radius
        self.rise = rise_p_residue
        self.omega = (2*np.pi)/residues_p_turn #Not hardcoding 1.74532925rad as engine should have more freedom

    def generate_helix(self,num_residues,centre_offset=(0.0,0.0,0.0)):

        x_start, y_start, z_start = centre_offset

        coordinates = []

        for i in range(num_residues):
            x = x_start + (self.r*np.cos(self.omega*i))
            y = y_start + (self.r*np.sin(self.omega*i))
            z = z_start + (self.rise*i)

            coordinates.append([x,y,z])

        return np.array(coordinates)

    
helix = ProteinBackboneGenerator().generate_helix(10)
print(helix)


       
        