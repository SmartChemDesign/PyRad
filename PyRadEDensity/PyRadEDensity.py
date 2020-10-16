import copy
import os

import matplotlib.pyplot as plt
import numpy as np
from sys import argv
from matplotlib.widgets import Slider, Button
from mpl_toolkits import mplot3d
from scipy import ndimage


class Cube:
    """
    Cube file class.
    For more broad tools check Cube-Toolz repo: https://github.com/funkymunkycool/Cube-Toolz
    """
    def __init__(self, _name=None):
        self.filename = ''


        #Two lines of text at the head of the file.
        self.comment_1 = 'none'
        self.comment_2 = 'none'

        #Defines the number of rows of molecular geometry data
        self.n_atoms = 0

        """
        This set of three fields defines the displacement vector from the geometric origin of the system 
        (0,0,0) to the reference point (x0,y0,z0) for the spanning vectors defined in {XAXIS}, {YAXIS}, and {ZAXIS}.
        """
        self.origin = np.array([0, 0, 0])

        """
        The first field on this line is an integer indicating the number of voxels NX present along the 
        X-axis of the volumetric region represented by the CUBE file. This value SHOULD always be positive; 
        whereas the input to the cubegen [Gau16] utility allows a negative value here as a flag for the units 
        of the axis dimensions, in a CUBE file distance units MUST always be in Bohrs, and thus the ‘units flag’
        function of a negative sign is superfluous. It is prudent to design applications to handle gracefully 
        a negative value here, however.

        The second through fourth values on this line are the components of the vector X⃗  defining the voxel X-axis. 
        They SHOULD all be non-negative; proper loading/interpretation/calculation behavior is not guaranteed if 
        negative values are supplied. As noted in the Gaussian documentation [Gau16], the voxel axes need neither 
        be orthogonal nor aligned with the geometry axes. However, many tools only support voxel axes that are 
        aligned with the geometry axes (and thus are also orthogonal). In this case, the first float value (Xx) 
        will be positive and the other two (Xy and Xz) will be identically zero.
        """
        self.n_x = 0
        self.n_y = 0
        self.n_z = 0

        self.x = 0
        self.y = 0
        self.z = 0

        # list of atom types
        self.atoms = []
        # list of atom positions
        self.atoms_xyz = []

        # data array
        self.data = []

        if _name is not None:
            try:
                self.load(_name)
            except IOError:
                print(f"Can't open {_name} file!")
                raise

    def load(self, _name):
        """
        Method to read cube file. Just needs the filename
        """
        with open(_name, 'r') as fin:
            self.filename = _name

            self.comment_1 = fin.readline()  # Save 1st comment
            self.comment_2 = fin.readline()  # Save 2nd comment

            _str = fin.readline().split()  # Number of Atoms and Origin
            self.n_atoms = int(_str[0])  # Number of Atoms
            self.origin = np.array([float(_str[1]), float(_str[2]), float(_str[3])])  # Position of Origin

            nVoxel = fin.readline().split()  # Number of Voxels
            self.n_x = int(nVoxel[0])
            self.x = np.array([float(nVoxel[1]), float(nVoxel[2]), float(nVoxel[3])])

            nVoxel = fin.readline().split()  #
            self.n_y = int(nVoxel[0])
            self.y = np.array([float(nVoxel[1]), float(nVoxel[2]), float(nVoxel[3])])

            nVoxel = fin.readline().split()  #
            self.n_z = int(nVoxel[0])
            self.z = np.array([float(nVoxel[1]), float(nVoxel[2]), float(nVoxel[3])])

            self.atoms = []
            self.atoms_xyz = []
            for atom in range(self.n_atoms):
                line = fin.readline().split()
                self.atoms.append(line[0])
                self.atoms_xyz.append(list(map(float, [line[2], line[3], line[4]])))

            self.data = np.zeros((self.n_x, self.n_y, self.n_z))

            i = int(0)
            for s in fin:
                for v in s.split():
                    self.data[int(i / (self.n_y * self.n_z)), int((i / self.n_z) % self.n_y),
                              int(i % self.n_z)] = float(v)
                    i += 1

        return None

    def save(self, _name):
        """
        Write out a Gaussian Cube file
        """
        try:
            with open(_name, 'w+') as fout:
                fout.write(".cube file generated from prt_esolv.py\n")
                fout.write(f"{_name}\n")

                fout.write(
                    f"{int(self.n_atoms)} {float(self.origin[0])} {float(self.origin[1])} {float(self.origin[2])}\n")

                fout.write(f"{int(self.n_x)} {float(self.x[0])} {float(self.x[1])} {float(self.x[2])}\n")
                fout.write(f"{int(self.n_y)} {float(self.y[0])} {float(self.y[1])} {float(self.y[2])}\n")
                fout.write(f"{int(self.n_z)} {float(self.z[0])} {float(self.z[1])} {float(self.z[2])}\n")

                for atom, xyz in zip(self.atoms, self.atoms_xyz):
                    fout.write(f"{atom} 0 {xyz[0]} {xyz[1]} {xyz[2]}\n")

                for ix in range(self.n_x):
                    for iy in range(self.n_y):
                        for iz in range(self.n_z):
                            fout.write(f"{self.data[ix][iy][iz]}")
                            if iz % 6 == 5:
                                fout.write('\n')
                        fout.write("\n")
        except IOError:
            print(f"Can't create {_name} file!!!")
            raise

        return None

    def info(self):
        """
        Writes out cube info
        """
        print(f"filename: {self.filename}")
        print(f"comments: \n{self.comment_1}{self.comment_2}")
        print(f"origin: {self.origin[0]}, {self.origin[1]}, {self.origin[2]}")
        print(f"atoms count: {self.n_atoms}")
        print(f"voxels count: {self.n_x}, {self.n_y}, {self.n_z}")
        print(f"voxel x-axis: {self.x[0]}, {self.x[1]}, {self.x[2]}")
        print(f"voxel y-axis: {self.y[0]}, {self.y[1]}, {self.y[2]}")
        print(f"voxel z-axis: {self.z[0]}, {self.z[1]}, {self.z[2]}")

    @property
    def xyz(self):
        return {"el": self.atoms, "xyz": self.atoms_xyz}

    def norm(self):
        """
        Shifts data origin to zero
        :return: None
        """
        old_origin = np.array(self.origin)
        self.origin = [0, 0, 0]
        old_origin[0] = old_origin[0] / self.x[0]
        old_origin[1] = old_origin[1] / self.y[1]
        old_origin[2] = old_origin[2] / self.z[2]
        self.data = ndimage.shift(self.data, -old_origin, mode='wrap')

    def density_maxima(self, samplesize=5, thresh_mod=0):
        """
        Gets all maxima from electronic density
        :param samplesize: size of summing volume
        :param thresh_mod: sensitivity of method
        :return: 3d array of maxima coordinates and 1d array of its value
        """
        filtered = ndimage.maximum_filter(self.data, size=(samplesize, samplesize, samplesize), mode="wrap")

        threshold = filtered.mean() + thresh_mod
        print(f"actual threhold value: {threshold:.2f}")
        labels, num_labels = ndimage.label(filtered > threshold)

        # Coordinates of maxima
        pos = np.array(ndimage.measurements.center_of_mass(np.asarray(self.data), labels=labels,
                                                           index=np.arange(1, num_labels + 1)))

        # Values of maxima
        val = np.array(ndimage.measurements.maximum(self.data, labels=labels, index=np.arange(1, num_labels + 1)))

        pos[:, 0] *= iCube.x[0]
        pos[:, 1] *= iCube.y[1]
        pos[:, 2] *= iCube.z[2]

        return pos, val


# MAIN ROUTINE #########################################################################################################
# check arguments line
if len(argv) == 1:
    print("Argument: .cube file")
    exit()

# process cube file
iCube = Cube(argv[1])
iCube.norm()
mpos, mval = iCube.density_maxima()

# plot init routine
fig = plt.figure()
ax = fig.add_subplot(111, projection=mplot3d.Axes3D.name)
plt.subplots_adjust(left=0.25, bottom=0.25)
ax.margins(x=0)

# setting up sliders
axcolor = 'lightgoldenrodyellow'
ax_ss = plt.axes([0.25, 0.19, 0.65, 0.03], facecolor=axcolor)
ax_th = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_corr = plt.axes([0.25, 0.11, 0.65, 0.03], facecolor=axcolor)
ax_save = plt.axes([0.25, 0.07, 0.65, 0.03])

S_ss = Slider(ax_ss, "sample size", 1, 10.0, valinit=5, valstep=1)
S_th = Slider(ax_th, "threshold", -3.5, 3.5, valinit=0, valstep=0.1)
S_corr = Slider(ax_corr, "correlation", 0, 8, valinit=0, valstep=0.5)
B_save = Button(ax_save, "Save anomalies")


def update(_):
    global o_type
    global img
    global cbar
    global mpos
    global mval

    cell_size = [iCube.x[0] * iCube.n_x, iCube.y[1] * iCube.n_y, iCube.z[2] * iCube.n_z]

    mpos, mval = iCube.density_maxima(int(S_ss.val), S_th.val)

    exclude = list()
    __maxima = np.asarray([0, 0, 0])

    if S_corr.val != 0:
        for i, _maxima in zip(range(len(mpos)), mpos):
            for pbc_x in range(-1, 1, 1):
                for pbc_y in range(-1, 1, 1):
                    for pbc_z in range(-1, 1, 1):
                        __maxima[0] = _maxima[0] + pbc_x * cell_size[0]
                        __maxima[1] = _maxima[1] + pbc_y * cell_size[1]
                        __maxima[2] = _maxima[2] + pbc_z * cell_size[2]
                        for _atomcoo in iCube.atoms_xyz:
                            if np.linalg.norm(_atomcoo - __maxima) < S_corr.val:
                                exclude.append(i)
                                break

        mpos = np.delete(mpos, exclude, axis=0)
        mval = np.delete(mval, exclude, axis=0)

        print(f"reduced maxima count: {len(mpos)}")
    else:
        print(f"maxima count: {len(mpos)}")

    img.remove()
    img = ax.scatter(mpos[:, 0], mpos[:, 1], mpos[:, 2], c=mval, cmap=plt.viridis())
    cbar.set_clim(vmin=mval.min(), vmax=mval.max())
    cbar.draw_all()
    fig.canvas.draw_idle()


# set update function for slider change
S_ss.on_changed(update)
S_th.on_changed(update)
S_corr.on_changed(update)


def save_anomalies(_):
    global mpos
    global mval
    global iCube

    anoCube = copy.copy(iCube)

    anoCube.n_atoms = len(mpos)
    anoCube.atoms = list(map(int, mval * 10))
    anoCube.atoms_xyz = list(mpos)
    
    _, tail = os.path.split(anoCube.filename)
    anoCube.save(tail.split(".")[0] + "_anomaly.cube")
    print("anomaly file saved!")


# set click action
B_save.on_clicked(save_anomalies)

# draw atoms
h_type = []
o_type = []
a_type = []

# sort atoms by type
for coo, t in zip(iCube.atoms_xyz, iCube.atoms):
    t = int(t)
    if t == 1:
        h_type.append(coo)
    elif t == 2:
        a_type.append(coo)
    else:
        o_type.append(coo)

# draw atoms
ax.plot(np.asarray(h_type)[:, 0], np.asarray(h_type)[:, 1], np.asarray(h_type)[:, 2], marker='$H$', color='black',
        linestyle="None")
ax.plot(np.asarray(a_type)[:, 0], np.asarray(a_type)[:, 1], np.asarray(a_type)[:, 2], marker='$a$', color='black',
        linestyle="None")
ax.plot(np.asarray(o_type)[:, 0], np.asarray(o_type)[:, 1], np.asarray(o_type)[:, 2], marker='$O$', color='black',
        linestyle="None")

# draw electron density maxima
img = ax.scatter(mpos[:, 0], mpos[:, 1], mpos[:, 2], c=mval, cmap=plt.viridis())
cbar = plt.colorbar(img, cax=fig.add_axes([0.1, 0.1, 0.03, 0.8]))

plt.show()
