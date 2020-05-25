# A little package for computing Voronoi diagrams within bounded regions using scipy
#
# A lot of this is from the help of a Stack Overflow 
# post I can no longer find with slight modifications
#
# How to use
# voronoi() to generate a diagram for the given input points in a bounded region
# generate_CVD() to generate a centroidal voronoi diagram within a bounded box with x iterations of loyd's algorithm
# plot_vornoi_diagram() plot a bounded voronoi diagram

import matplotlib.pyplot as pl
import numpy as np
import scipy as sp
import scipy.spatial
import sys
import copy

eps = sys.float_info.epsilon

# returns a new np.array of towers that within the bounding_box
def in_box(towers, bounding_box):
    return np.logical_and(np.logical_and(bounding_box[0] <= towers[:, 0],
                                         towers[:, 0] <= bounding_box[1]),
                          np.logical_and(bounding_box[2] <= towers[:, 1],
                                         towers[:, 1] <= bounding_box[3]))


#generates a bounded vornoi diagram with finite regions
def voronoi(towers, bounding_box):
    # Select towers inside the bounding box
    i = in_box(towers, bounding_box)
    # Mirror points
    points_center = towers[i, :]
    points_left = np.copy(points_center)
    points_left[:, 0] = bounding_box[0] - (points_left[:, 0] - bounding_box[0])
    points_right = np.copy(points_center)
    points_right[:, 0] = bounding_box[1] + (bounding_box[1] - points_right[:, 0])
    points_down = np.copy(points_center)
    points_down[:, 1] = bounding_box[2] - (points_down[:, 1] - bounding_box[2])
    points_up = np.copy(points_center)
    points_up[:, 1] = bounding_box[3] + (bounding_box[3] - points_up[:, 1])
    points = np.append(points_center,
                       np.append(np.append(points_left,
                                           points_right,
                                           axis=0),
                                 np.append(points_down,
                                           points_up,
                                           axis=0),
                                 axis=0),
                       axis=0)
    # Compute Voronoi
    vor = sp.spatial.Voronoi(points)
    # Filter regions
    regions = []
    for region in vor.regions:
        flag = True
        for index in region:
            if index == -1:
                flag = False
                break
            else:
                x = vor.vertices[index, 0]
                y = vor.vertices[index, 1]
                if not(bounding_box[0] - eps <= x and x <= bounding_box[1] + eps and
                       bounding_box[2] - eps <= y and y <= bounding_box[3] + eps):
                    flag = False
                    break
        if region != [] and flag:
            regions.append(region)
    vor.filtered_points = points_center # creates a new attibute for points that form the diagram within the region
    vor.filtered_regions = regions
    return vor


# finds the centroid of a region
def centroid_region(vertices):
    # Polygon's signed area
    A = 0
    # Centroid's x
    C_x = 0
    # Centroid's y
    C_y = 0
    for i in range(0, len(vertices) - 1):
        s = (vertices[i, 0] * vertices[i + 1, 1] - vertices[i + 1, 0] * vertices[i, 1])
        A = A + s
        C_x = C_x + (vertices[i, 0] + vertices[i + 1, 0]) * s
        C_y = C_y + (vertices[i, 1] + vertices[i + 1, 1]) * s
    A = 0.5 * A
    C_x = (1.0 / (6.0 * A)) * C_x
    C_y = (1.0 / (6.0 * A)) * C_y
    return np.array([[C_x, C_y]])


# performs x iterations of loyd's algorithm to calculate a centroidal vornoi diagram
def generate_CVD(points, iterations, bounding_box):
    p = copy.copy(points)
    for i in range(iterations):
        # print("there are " + str(len(p)) + " towers in the vornoi diagram right now")
        
        vor = voronoi(p, bounding_box)
        centroids = []
        for region in vor.filtered_regions:
            vertices = vor.vertices[region + [region[0]], :]
            centroid = centroid_region(vertices)
            centroids.append(list(centroid[0, :]))

        p = np.array(centroids)

    return voronoi(p, bounding_box)
        

# plots vornoi data on a pyplot
def plot_vornoi_diagram(vor, bounding_box):
    # Initializes pyplot stuff
    fig = pl.figure()
    ax = fig.gca()

    # Plot initial points
    ax.plot(vor.filtered_points[:, 0], vor.filtered_points[:, 1], 'b.')

    # Plot ridges points
    for region in vor.filtered_regions:
        vertices = vor.vertices[region, :]
        ax.plot(vertices[:, 0], vertices[:, 1], 'go')

    # Plot ridges
    for region in vor.filtered_regions:
        vertices = vor.vertices[region + [region[0]], :]
        ax.plot(vertices[:, 0], vertices[:, 1], 'k-')

    # stores references to numbers for setting axes limits
    margin_percent = .1
    width = bounding_box[1]-bounding_box[0]
    height = bounding_box[3]-bounding_box[2]

    ax.set_xlim([bounding_box[0]-width*margin_percent, bounding_box[1]+width*margin_percent])
    ax.set_ylim([bounding_box[2]-height*margin_percent, bounding_box[3]+height*margin_percent])

    pl.show()