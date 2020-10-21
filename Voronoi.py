# A little package for computing Voronoi diagrams within bounded regions using scipy
#
# A lot of this is from the help of a Flabetvvibes and Energya on StackOverflow
# 
#
# How to use
# bounded_voronoi() to generate a diagram for the given input points in a bounded region
# bounded_CVD() to generate a centroidal voronoi diagram within a bounded box with x iterations of loyd's algorithm
# plot_vornoi_diagram() plot a bounded voronoi diagram

import matplotlib.pyplot as pl
import matplotlib.path as pyPath
from shapely.geometry.polygon import Polygon, Point
import shapely
import numpy as np
import scipy as sp
import scipy.spatial
import sys
import copy

eps = sys.float_info.epsilon

# Returns a new np.array of towers that within the bounding_box
def in_box(towers, bounding_box):
    return np.logical_and(np.logical_and(bounding_box[0] <= towers[:, 0],
                                         towers[:, 0] <= bounding_box[1]),
                          np.logical_and(bounding_box[2] <= towers[:, 1],
                                         towers[:, 1] <= bounding_box[3]))


# Generates a bounded vornoi diagram with finite regions
def bounded_voronoi(towers, bounding_box):
    # Select towers inside the bounding box
    i = in_box(towers, bounding_box)

    # Mirror points left, right, above, and under to provide finite regions for the edge regions of the bounding box
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

    vor.filtered_points = points_center # creates a new attibute for points that form the diagram within the region
    vor.filtered_regions = np.array(vor.regions)[vor.point_region[:vor.npoints//5]] # grabs the first fifth of the regions, which are the original regions

    return vor


# Finds the centroid of a region. First and last point should be the same.
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


# Performs x iterations of loyd's algorithm to calculate a centroidal vornoi diagram
def bounded_CVD(points, iterations, bounding_box):
    p = copy.copy(points)

    for i in range(iterations):
        vor = bounded_voronoi(p, bounding_box)
        centroids = []

        for region in vor.filtered_regions:
            vertices = vor.vertices[region + [region[0]], :] # grabs vertices for the region and adds a duplicate of the first one to the end
            centroid = centroid_region(vertices)
            centroids.append(list(centroid[0, :]))

        p = np.array(centroids)

    return bounded_voronoi(p, bounding_box)


# returns a cummulative summation matrix of given 2d array
def generate_cummulative_summation_matrix(weights):
    width, height = weights.shape
    cumsum = np.zeros(weights.shape)

    for x in range (1,width):
        for y in range(height):
            cumsum[x, y] = cumsum[x-1, y] + weights[x, y]

    return cumsum


def get_Bounding_Box(vertices):
    bb = [vertices[0][0], 0, vertices[0][1], 0] # [xmin, xmax, ymin, ymax]

    for vertex in vertices:
        if(vertex[0] < bb[0]):
            bb[0] = vertex[0]
        if(vertex[0] > bb[1]):
            bb[1] = vertex[0]
        if(vertex[1] < bb[2]):
            bb[2] = vertex[1]
        if(vertex[1] > bb[3]):
            bb[3] = vertex[1]

        bb[0] = np.floor(bb[0])-1
        bb[2] = np.floor(bb[2])-1
        bb[1] = np.ceil(bb[1])+1
        bb[3] = np.ceil(bb[3])+1
    
    return bb

# returns coordinates of the weighted centroid using Grant Trebbin's techinique
def calc_weighted_centroid(P, Q, vertices):
    width, height = P.shape
    vertices = [ [i[0]-1, i[1]-1] for i in vertices] # adjust x values for vertices to compensate for zero padding in P/Q

    px = []
    polygon = Polygon(vertices)
    bounding_box = get_Bounding_Box(vertices)

    for y in range(int(bounding_box[2]), int(bounding_box[3])):
        row = []
        for x in range(int(bounding_box[0]), int(bounding_box[1])):
            if polygon.contains(Point(x, y)) or polygon.touches(Point(x, y)):
                row.append((x, y))
        if(row != []):
            px.append(row)


    denom = 0
    xNum = 0
    yNum = 0

    for row in px:
        first = row[0]
        last = row[-1]

        denom += P[last] - P[first]
        xNum += (last[0]*P[last]-first[0]*P[first]) - (Q[last]-Q[first])
        yNum += (width-last[1])*(P[last]-P[first])

    if(denom == 0 or xNum == 0 or yNum == 0): # todo denom, xNum, and yNum become zero when point density is higher
        print("trying to divide by zero")
        print("bounding box for cell: " + str(bounding_box))
        # print("pixels in cell: " + str(px))

    cx = xNum/denom
    cy = yNum/denom
    # print("denom = " + str(denom) + " xNum = " + str(xNum) + " yNum = " + str(yNum) + " centroid = (" + str(cx) + ", " + str(cy) + ")")

    return np.array([[cx, cy]])


# Performs x iterations of loyd's algorithm to calculate a weighted centroidal vornoi diagram
def bounded_weighted_CVD(points, iterations, bounding_box, weights):
    print("Voronoi.py:: creating bounded_weighted_CVD")
    pt = copy.copy(points)


    weights = np.pad(weights, ((1, 0), (0, 0)), mode='constant') # zero pads to prep for cummulative summation
    print("Voronoi.py:: making first Summation Matrix P")
    P = generate_cummulative_summation_matrix(weights)
    print("Voronoi.py:: making second Summation Matrix Q")
    Q = generate_cummulative_summation_matrix(P)

    for i in range(iterations):
        print(len([x for x in pt if x[0] < bounding_box[0] or x[0] > bounding_box[1] or x[1] < bounding_box[2] or x[1] > bounding_box[3]]))
        vor = bounded_voronoi(pt, bounding_box)
        print("Voronoi.py:: " + str(i+1) + " iteration of " + str(iterations) + " of loyd's alogithm:: " + str(len(vor.filtered_regions)) + " regions in diagram")
        centroids = []

        for region in vor.filtered_regions:
            vertices = vor.vertices[region] # grabs vertices for the region and adds a duplicate of the first one to the end
            centroid = calc_weighted_centroid(P, Q, vertices)
            centroids.append(list(centroid[0, :]))

        pt = np.array(centroids)

    return bounded_voronoi(pt, bounding_box)
        
# returns a pyplot of given voronoi data
def plot_vornoi_diagram(vor, bounding_box, show_figure):
    print("Voronoi.py:: Plotting Voronoi Diagram")
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

    if show_figure:
        pl.show()
    return fig