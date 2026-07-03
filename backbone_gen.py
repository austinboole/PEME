import numpy as np
import matplotlib.pyplot as plt

class ProteinBackboneGenerator:
    def __init__(self,radius=2.3,rise_p_residue=1.5,residues_p_turn=3.6):
        self.r = radius
        self.rise = rise_p_residue
        self.omega = (2*np.pi)/residues_p_turn #Not hardcoding 1.74532925rad as engine should have more freedom

    def generate_helix(self,num_residues,centre_offset=(0.0,0.0,0.0),direction=1):
        x_start, y_start, z_start = centre_offset
        coordinates = []


        helix_height = self.rise*(num_residues-1) #Calculating the maximum height of the helix (so that helices can be drawn downward)

        for i in range(num_residues):
            x = x_start + (self.r*np.cos(self.omega*i))
            y = y_start + (self.r*np.sin(self.omega*i))


            z = z_start + helix_height*((1-direction)/2) + direction*self.rise*i #Z axis increases with residues when direction=1, but decreases when direction=-1

            coordinates.append([x,y,z])

        return np.array(coordinates)
    
    def generate_bundle(self, num_residues):
        #Calculate symmetrical grid
        bundle_helices = []
        bundle_radius = 5.0 #distance from origin to centre of any helix

        for k in range(4):
            theta = (np.pi/4) + ((k*np.pi)/2) #Calculating angle of 45,135,225,315 for the placement of each helix

            offset = (bundle_radius*np.cos(theta), bundle_radius*np.sin(theta), 0.0)
            quadrant_direction = (-1)**k

            helix = self.generate_helix(num_residues,centre_offset=offset,direction=quadrant_direction)
            bundle_helices.append(helix)
        return bundle_helices


bundle = ProteinBackboneGenerator().generate_bundle(20)
print(bundle)

#3D map visualisation

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')
for index, helix in enumerate(bundle):
    x_coords = helix[:, 0]
    y_coords = helix[:, 1]
    z_coords = helix[:, 2]

    ax.plot(x_coords,y_coords,z_coords)
ax.set_box_aspect([1,1,1])
plt.show()