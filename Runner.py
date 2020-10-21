from PIL import Image
from Stipple import stipple
import random
import copy

im = Image.open('exampleR.png') # loads in an image
# im = Image.open('pictureSQ.jpg') # loads in an image
# im = Image.open('pictureRect.jpg') # loads in an image
stipple(im, True, 100, 20) # creates a stippled version of the image

