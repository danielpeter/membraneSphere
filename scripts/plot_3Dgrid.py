#!/usr/bin/env python
#
# plot the cross-correlation and L2-norm between reference and output seismograms
#
from __future__ import print_function

import sys
import os

import numpy as np

# matplotlib
import matplotlib.pyplot as plt
# 3D plotting
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# coloring
from matplotlib import cm
from matplotlib.colors import Normalize,LightSource


def plot_grid(Dvvert_file,Dvface_file,cell_file):
    """
    plots spherical grid
    """
    print("Spherical grid:")
    print(f"  vertices file : {Dvvert_file}")
    print(f"  faces    file : {Dvface_file}")
    if cell_file:
        print(f"  cell data file: {cell_file}")
    print("")

    # checks if files exists
    if not os.path.isfile(Dvvert_file):
        print("Please check if file exists: ",Dvvert_file)
        sys.exit(1)
    if not os.path.isfile(Dvface_file):
        print("Please check if file exists: ",Dvface_file)
        sys.exit(1)

    # check if cell data provided
    use_cell_data = False
    if cell_file:
        # checks if files exists
        if not os.path.isfile(cell_file):
            print("Please check if file exists: ",cell_file)
            sys.exit(1)

        # use data
        use_cell_data = True

    # gets file data
    # vertex coordinates
    # format: #x #y #z
    Dvvert = np.loadtxt(Dvvert_file, dtype=np.float64)

    # face vertex ids
    # format: #id1 #id2 .. #id6
    Dvface = np.loadtxt(Dvface_file, dtype=np.integer)

    # Adjust Dvface to be 0-based indexing
    Dvface = Dvface - 1

    print("vertices : ",len(Dvvert[:,0]))
    print("faces    : ",len(Dvface[:,0]))
    print("")

    # cell data
    if use_cell_data:
        # format: #val
        cellData = np.loadtxt(cell_file, dtype=np.float64)

        # cell data filename as label
        # Extract the basename (removes directory path)
        basename = os.path.basename(cell_file)  # Gets 'cellFractionAverage.dat'
        # Remove the extension
        cellDataName = os.path.splitext(basename)[0]  # Gets 'cellFractionAverage'

        print("cell data:")
        print(f"  length   : {len(cellData)}")
        print(f"  min/max  : {cellData.min()} / {cellData.max()}")
        print(f"  data name: {cellDataName}")
        print("")

        # check data length
        if len(cellData) != len(Dvface[:,0]):
            print(f"Error: mismatch in cell data length - should be same as number of faces {len(Dvface[:,0])}")
            print( "       Please check, exiting...")
            sys.exit(1)

        # Choose a colormap (you can select from matplotlib.cm)
        cmap = cm.viridis

        # Normalize the data to the range [0, 1] for the colormap
        norm = Normalize(vmin=cellData.min(), vmax=cellData.max())

    # labelling infos
    # Extract the basename from Dvvert file
    basename = os.path.basename(Dvvert_file)  # Gets 'Dvvert4.dat'
    # get subdivisions number from file name
    subdivisions = ''.join(filter(str.isdigit, basename))

    print("grid subdivisions : ",subdivisions)
    print("")

    # 3D plotting
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection='3d')

    # Prepare the vertices for the Poly3DCollection
    face_vertices = []
    face_colors = []
    for i, face in enumerate(Dvface):
        # Get the vertices for the current face
        #print(f"{i} face {face}")
        vertices = Dvvert[face]
        face_vertices.append(vertices)
        # coloring
        color = 'lightblue'
        # cell data coloring
        if use_cell_data:
            # Get the data value for the current face
            cvalue = cellData[i]
            # Map the data value to a color using the colormap and normalization
            color = cmap(norm(cvalue))
        face_colors.append(color)

    face_vertices = np.array(face_vertices)
    face_colors = np.array(face_colors)

    # Create the Poly3DCollection
    if use_cell_data:
        # light
        ls = LightSource(260, -80)  # azimuth from north, elevation from zero plane in degrees
        # polygon object
        poly = Poly3DCollection(face_vertices, facecolors=face_colors, edgecolor='black',linewidth=0.1,
                                shade=True,lightsource=ls)
        # set color limits on poly object
        poly.set_clim(cellData.min(),cellData.max())
    else:
        poly = Poly3DCollection(face_vertices, facecolors=face_colors, edgecolor='black',linewidth=0.1)

    # Add the collection to the axes
    ax.add_collection3d(poly)

    # Auto-scale the view to the data
    all_vertices = np.vstack(Dvvert[Dvface])
    min_x, max_x = np.min(all_vertices[:, 0]), np.max(all_vertices[:, 0])
    min_y, max_y = np.min(all_vertices[:, 1]), np.max(all_vertices[:, 1])
    min_z, max_z = np.min(all_vertices[:, 2]), np.max(all_vertices[:, 2])
    ax.set_xlim(min_x - 0.1, max_x + 0.1)
    ax.set_ylim(min_y - 0.1, max_y + 0.1)
    ax.set_zlim(min_z - 0.1, max_z + 0.1)
    ax.set_aspect('equalxy')

    # Add a colorbar to show the mapping of data to colors
    if use_cell_data:
        # sm - not working in 3D...
        ##sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        ##sm.set_array([])  # For older versions of matplotlib
        ##cbar = fig.colorbar(sm)
        # poly has wrong scale [0,1] since color values are normed
        cbar = fig.colorbar(poly, shrink=0.5, aspect=20, extend='both')
        # cell data filename as label
        cbar.set_label(f"{cellDataName}")

    # annotation
    plt.title(f"Spherical grid - subdivisions {subdivisions}", fontsize=16, fontweight='bold')
    # Add a subtle background style (optional)
    #plt.style.use('ggplot')
    # Adjust layout to prevent labels from overlapping
    #plt.tight_layout()

    # saves as JPEG file
    filename = "out.jpg"
    plt.savefig(filename)
    print("plotted to: out.jpg")
    print("")

    # show figure
    plt.show()


def usage():
    print("usage: ./plot_3Dgrid.py DvvertN.dat DvfaceN.dat [cell.dat]")
    print("  with")
    print("     DvvertN.dat - vertex file, e.g. griddata/Dvvert4.dat")
    print("     DvfaceN.dat - face file, e.g. griddata/Dvface4.dat")
    print("     cell.dat    - (optional) corresponding cell data, e.g. OUTPUT/cellFractionAverage.dat")
    sys.exit(1)

if __name__ == '__main__':
    # input
    cell_file = None
    if len(sys.argv) < 3: usage()

    # gets arguments
    Dvvert_file = sys.argv[1]
    Dvface_file = sys.argv[2]
    if len(sys.argv) >= 4: cell_file = sys.argv[3]

    plot_grid(Dvvert_file,Dvface_file,cell_file)

