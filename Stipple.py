from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
from Voronoi import bounded_voronoi, generate_CVD, plot_vornoi_diagram
import numpy as np
import random
import math
import copy
import time

# plots a stippled verison of the image
def stipple(original_image, importance_sampling, stipple_amount, iterations):
    # stores a gray scale copy of the image
    image = original_image.convert('LA')

    points = generate_initial_points(image, stipple_amount, importance_sampling)
    bounding_box = [0, image.size[0], 0, image.size[1]]

    # generates centroidal voronoi diagram
    vor = generate_CVD(points, iterations, bounding_box)

    # plots voronoi and displays image
    fig = plot_vornoi_diagram(vor, bounding_box, False)
    plt.imshow(original_image,zorder=0)
    plt.show()


# generates points randomly or by using importance sampling for the given image
def generate_initial_points(image, stipple_amount, importance_sampling):
    width, height = image.size
    pixels = image.load()
    points = []

    if(importance_sampling): # creates initial set of points for stippling randomly
        possible_points = []
        for x in range(width):
            for y in range(height):
                possible_points.append([x,y])
        
        for x in range(stipple_amount):
            r = np.random.randint(len(possible_points))
            points.append(possible_points[r])
            possible_points.pop(r)
    
    else: # creates initial set of points for stippling with "importance sampling" from Anton Lopyrev
        darkness_histogram = []
        for x in range(256):
            darkness_histogram.append([])

        for x in range(width): # puts pixels in histogram
            for y in range(height):
                darkness_histogram[pixels[x, y][0]].append([x, y])

        weighted_bins = []
        last_weight = 0

        for bin_num in range(len(darkness_histogram)-5): # creates weights for pixels
            weight = (255 - bin_num)*len(darkness_histogram[bin_num])
            weighted_bins.append(last_weight + weight)
            last_weight += weight

        for x in range(stipple_amount): # generates a random number and finds the bin and pixel it corresponds to
            r = np.random.randint(0, last_weight)
            random_bin = find_weighted_bin(r, weighted_bins)
            r = np.random.randint(0, len(darkness_histogram[random_bin]))
            points.append(darkness_histogram[random_bin][r])

    return np.array(points)


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