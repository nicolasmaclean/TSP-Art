from PIL import Image
from Stipple import stipple
import random
import copy

im = Image.open('picture.jpg') # loads in an image
stipple(im, True, 100, 5) # creates a stippled version of the image

