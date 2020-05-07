import cv2
import numpy as np
from matplotlib import pyplot as plt
from scipy import ndimage
from skimage import measure, color, io
from skimage.segmentation import clear_border

def segment(original_image, image, pixels_to_um, cell_or_nucleus_label, file_path):
    #Image converted to binary
    ret, thresh = cv2.threshold(image, 0, 255, cv2.THRESH_BINARY+cv2.THRESH_OTSU)
    
    #Removes small noise
    kernel_1 = np.ones((3,3), np.uint8)
    opening = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel_1, iterations = 1)
    
    #Removes the cells/nuclei on edges of the image and identifies the background as not part of cell or nuclei
    opening = clear_border(opening)
    sure_bg = cv2.dilate(opening, kernel_1, iterations = 10)
    
    #Finding the required area the difference between the background and cell/nuclei boundary
    dist_transform = cv2.distanceTransform(opening, cv2.DIST_L2, 3)
    ret2, sure_fg = cv2.threshold(dist_transform, 0.01*dist_transform.max(), 255, 0)
    
    #Unknown and ambiguous background is removed
    sure_fg = np.uint8(sure_fg)
    unknown = cv2.subtract(sure_bg, sure_fg)
    
    #markers created to label regions and required regions are labelled as positive numbers
    ret3, markers = cv2.connectedComponents(sure_fg)
    markers += 10
    
    #All unrequired markers are set to 0
    markers[unknown == 255] = 0
    
    #Required markers are added on top of original image to highlight the difference between nuclei and cells
    markers = cv2.watershed(original_image, markers)
    original_image[markers == -1] = [0,255,255]
    
    #Maintaing original RGB colour scheme after watershed function and colours in yellow and blue boundaries
    image2 = color.label2rgb(markers, bg_label = 0)
    regions = measure.regionprops(markers, intensity_image=image)
    
    #Final edited image is saved
    cv2.imwrite(file_path + ".jpg", original_image)
    
    return len(regions)

def main():
    file_path_green = "/Users/jawahar/Desktop/Programming/BioPython/Viability images/LPA/Image 23_Maximum intensity projection.tif"
    file_path_red = "/Users/jawahar/Desktop/Programming/BioPython/Viability images/LPA/Image 23_Maximum intensity projection.tif"
    image = cv2.imread("/Users/jawahar/Desktop/Programming/BioPython/Viability images/LPA/Image 23_Maximum intensity projection.tif")
    print(image)
    green, red = image[:, :, 1], image[:, :, 2]
    print(green)
    print(red)
    pixels_to_um = 1.2044

    #Changes the contrast and brightness of the YAP channel image to identify the cells
    new_green = np.zeros(green.shape, green.dtype)
    contrast = 3.0 #Contrast: from 1-3
    brightness = 100 #Brightness: from 0-100
    for pixel in range(green.shape[0]):
        new_green[pixel] = np.clip(contrast*green[pixel] + brightness, 0, 255)
    
    print(segment(image, new_green, pixels_to_um, "Green", file_path_green))
    print(segment(image, red, pixels_to_um, "Red", file_path_red))

main()
