from PIL import Image, ImageDraw
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy as np
import random
import math
import copy
import time


# plots a stippled verison of the image
def stipple(original_image, random_input, stipple_amount):
    im_copy = original_image.convert('LA') # stores a gray scale copy of the image
    width, height = im_copy.size # stores image dimensions
    im_copy_pixels = im_copy.load() # 2d array that contains references to the pixels of im

    points = []
    probability = .01

    if(random_input): # creates initial set of points for stippling randomly
        possible_points = []
        for x in range(width):
            for y in range(height):
                possible_points.append([x,y])
        
        for x in range(stipple_amount):
            r = np.random.randint(len(possible_points))
            points.append(possible_points[r])
            possible_points.pop(r)
    
    else: # creates initial set of points for stippling with "importance sampling" from Anton Lopyrev
        # darkness_histogram = [list().copy()]*256 # a histogram of pixels based by darkness value from grayscale image
        darkness_histogram = []
        for x in range(256):
            darkness_histogram.append([])

        for x in range(width): # puts pixels in histogram
            for y in range(height):
                darkness_histogram[im_copy_pixels[x, y][0]].append([x, y])

        weighted_bins = []
        last_weight = 0

        for bin_num in range(len(darkness_histogram)-5): # creates weights for pixels
            weight = (255 - bin_num)*len(darkness_histogram[bin_num])
            weighted_bins.append(last_weight + weight)
            last_weight += weight

        for x in range(stipple_amount):
            r = np.random.randint(0, last_weight)
            random_bin = find_weighted_bin(r, weighted_bins)
            r = np.random.randint(0, len(darkness_histogram[random_bin]))
            points.append(darkness_histogram[random_bin][r])
            
    # adds points just outside image to create defined regions for edge cells
    hack_points = [[-100, -100],
                   [-100, height+100],
                   [width+100, -100],
                   [width+100, height+100],

                   [-150, -150],
                   [-150, height+150],
                   [width+150, -150],
                   [width+150, height+150]]
    points += hack_points
            
    scatter_points(points, (width, height), original_image, hack_points)


# performs a binary search to find the corresponding weighted bin
def find_weighted_bin(num, weighted_bins):
    ls = weighted_bins.copy()
    ls.insert(0, -1)

    start = 1
    end = len(ls) - 1

    while start <= end:
        mid = int((start + end) / 2)
        if num == ls[mid] or (num < ls[mid] and num > ls[mid-1]):
            return mid
        elif num < ls[mid]:
            end = mid - 1
        elif num > ls[mid]:
            start = mid + 1
    return -1


    # for x in range(width): # generates stippling dots proportional to the darkness of each pixel
    #     for y in range(height):
    #         if(random.random()*255 < (255 - im_copy_pixels[x, y][0])*probability):
    #             points.append([x, y])
                

def scatter_points(input_points, image_size, img, hack_points):
    points = np.array(input_points) # uses input points

    plt.style.use('seaborn-dark') # sets up figure
    fig = plt.figure(1)
    axes = plt.gca()
    axes.set_xlim([0, image_size[0]]) # limits axes to the size of the input image
    axes.set_ylim([0, image_size[1]])

    plt.scatter(points[:, 0], points[:, 1]) # plots input points
    plt.imshow(img,zorder=0)

    start_time = time.time()

    vor = Voronoi(points) # generates vornoi points using scipy
    fig = voronoi_plot_2d(vor)

    axes = plt.gca()
    axes.set_xlim([0, image_size[0]])
    axes.set_ylim([0, image_size[1]])

    print(time.time()-start_time)
    start_time = time.time()

    vor = generate_cvd_points(vor, 5, hack_points)

    # fig = voronoi_plot_2d(vor) # plots vornoi points

    axes = plt.gca()
    axes.set_xlim([0, image_size[0]])
    axes.set_ylim([0, image_size[1]])

    print(time.time()-start_time)

    plt.show()

# generates points for a centroidal vornoi diagram 
def generate_cvd_points(vornoi, iterations, hack_points):
    vor = copy.copy(vornoi)

    for it in range(iterations):
        centroids = []

        for region in vor.regions:
            if not -1 in region and len(region) > 2:
                region_vertices = []

                for vertex in region:
                    region_vertices.append(vor.vertices[vertex])

                centroids.append(find_centroid(region_vertices))

        centroids += hack_points
        centroids = np.array(centroids) # uses input points
        print(len(centroids))

        vor = Voronoi(centroids) # generates vornoi points using scipy
        fig = voronoi_plot_2d(vor) # plots vornoi points

    return vor


# finds the centroid of a cell. Thanks to Emile Cormier on Stack Overflow and wikipedia
def find_centroid(vertices):
    centroid = [0, 0]
    signedArea = 0.0    # signed area
    x0 = 0.0            # Current vertex X
    y0 = 0.0            # Current vertex Y
    x1 = 0.0            # Next vertex X
    y1 = 0.0            # Next vertex Y
    a = 0.0             # Partial signed area

    # // For all vertices except last
    for i in range(len(vertices)-1):
        x0 = vertices[i][0] # todo fix .x and .y access
        y0 = vertices[i][1]
        x1 = vertices[i+1][0]
        y1 = vertices[i+1][1]
        a = x0*y1 - x1*y0
        signedArea += a
        centroid[0] += (x0 + x1)*a
        centroid[1] += (y0 + y1)*a

    # // Do last vertex separately to avoid performing an expensive
    # // modulus operation in each iteration.
    i = len(vertices)-1
    x0 = vertices[i][0] # todo also fix .x and .y access
    y0 = vertices[i][1]
    x1 = vertices[0][0]
    y1 = vertices[0][1]
    a = x0*y1 - x1*y0
    signedArea += a
    centroid[0] += (x0 + x1)*a
    centroid[1] += (y0 + y1)*a

    signedArea *= 0.5
    centroid[0] /= (6.0*signedArea)
    centroid[1] /= (6.0*signedArea)

    return centroid