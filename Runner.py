from PIL import Image
from Stipple import stipple
import random
import copy

im = Image.open('picture.jpg') # loads in an image
stipple(im, False, 100) # creates a stippled version of the image