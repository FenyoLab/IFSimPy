# IFSimPy
This repository contains python scripts implementing the calculation of the Interaction Factor (IF) for assessing the degree of interaction of molecules in fluorescent microscopy images.

----------------------------------------------------------------------------------------------------------------------------
interaction_factor.py

The class interaction_factor encapsulates all of the main functionality, including calculating the IF given an image and also producting simulations at a given IF level using the thresholded objects of the image.

Usage:

def __init__(self, rgb_mask_image, roi_mask=[], intensity_image=[])
The class is initialized by providing an rgb mask image, an optional ROI mask image and an optional intensity image.  The rgb mask image should be of type bool and shape n x m x 3.  Connected True values in the mask image represent the objects for the corresponding color channel.  This image can be produced by applying a threshold to each color channel separately and combining the corresponding mask images (image2mask.py provides this functionality).  This image can be produced by alternative means as well, such as a point based clustering method, typically from dSTORM/PALM data.  The ROI mask is of type bool and shape n x m, and outlines the area of the image where objects are expected to reside.  Simulations will be drawn in this ROI area.  The intensity image n x m x 3, is the original intensity image before thresholding. If given, the function providing simulations at a given IF level will output intensity images as well as mask images for the simulations.  (The pixels will only have an intensity value inside the segmented area of the objects.)  The following class variables are set or initialized to default values:
    ----------------------------------------------
    Image parameters:
    -----------------
    labeling_structure : used in labeling objects, default is 8-connectivity
    cluster_masks      : array, mask for each color channel
    num_clusters       : array, number of clusters in each color channel
    cluster_props      : array, cluster properties in each color channel
    roi_mask           : the roi, if using
    use_roi            : True if input roi_mask is not empty
    area_roi           : area of the ROI or entire image area if use_roi==False
    roi_min/max_row/col: min/max row/col of the bbox of the ROI 
                         (or of entire image), for internal use
    intensity_image    : intensity image, if using
    
    IF parameters:
    --------------
    nonref_color (0)  : nonreference color channel to use for simulations
    ref_color (1)     : reference color channel to use for simulations
    move_nonref (True): True/False whether to move the nonref color when making 
                        simulations
    
    Parameters used in calculate_IF:
    --------------------------------
    plot_IF_curve (False)     : saves a plot of the IF curve as given by the IF 
                                equation
    save_random_images (False): if True, random images produced to calculate IF 
                                will be saved
    num_random_sims (50)      : number of random sims to use to calcuated IF 
        
    Parameters used in simulate_IF:
    -------------------------------
    save_image (True)            : if True, simulated image will be saved
    save_colors (['red','green']): the colors to use if saving simulated image 
                                   [non-ref, ref], choices: 'red', 'green', 'blue', 'magenta' 
                                   (magenta should only be paired with green)
    image_file_prefix ('')       : prefix to use in file name of simulated image
    use_ellipses (False)         : if True, ellipses will be produced instead of 
                                   original clusters; the size/shape of ellipses 
                                   are taken from cluster sizes/shapes
    nonref_size (1)              : when use_ellipses==True, a number to multiply 
                                   by the size of the major/minor axis of the ellipses 
                                   of the nonreference color
    ref_size (1)                 : same as above for reference color
    nonref_number (None)         : if not None, the number of nonref clusters to 
                                   use in the simulation; this number of clusters 
                                   is drawn randomly from the cluster (or ellipse) 
                                   list, if this is set to None, the original nonref 
                                   cluser list from the image is used 
    ref_number (None)            : same as above for reference color
        
    Directory for saving output of either function:
    -----------------------------------------------
    save_dir (''): if is empty string, output not saved. 

def calculate_IF(self)
This function is called to calculate the IF for the image, for the non-reference/reference color channels specified (see parameters above).  The return value is (calcIF, p_val, orig_meas), where calcIF = the IF for the image, p_val = the fraction of random simulations where the %-overlap was >= %-overlap of the input image.  orig_meas = %-overlap for the input image.  Prior to running this function, change any class variables from their default values, as desired.  (see list above)

def simulate_IF(self, IF=0)
Produces simulated image at the given IF (IF=0 is random simulations), for the non-reference/reference color channels specified (see parameters above).  Set up the directory for the images in the class variable save_dir.  If desired, a prefix can be appended to the file name (image_file_prefix).  Prior to running this function, change any class variables from their default values, as desired (see list above).  Calculates the following measurements for the simulated image:
    sim_num_ref_ov_clusters: number of reference overlapping clusters in the simulation
    sim_overlap_mask       : overlap mask for ref+nonref clusters in the simulation
    sim_num_ov_clusters    : number of overlaps in overlap mask (simulation)
    sim_area_ov_clusters   : sum of total area of clusters in overlaop mask (simulation)
    simulated_images       : array, simulated image for each channel, includes intensities, 
                             if intensity image is not empty (otherwise values 
                             are either 0 or 255)
    simulated_masks        : array, simulated image for each channel,  boolean values
    ref_cluster_ov         : dict, T/F whether ref cluster with given label overlaps
                             a nonref cluster (label is from cluster_props[ref_color])

def calc_orig_image_measurements(self)
Calculates the following measurements for the input image:
    orig_area_clusters      : array, sum of total area of clusters, each channel
    orig_overlap_mask       : overlap mask for ref+nonref colors 
    orig_num_ov_clusters    : number of overlaps in overlap mask
    orig_area_ov_clusters   : sum of total area of clusters in overlap mask
    orig_num_ref_ov_clusters: number of reference clusters that overlap non-ref clusters

----------------------------------------------------------------------------------------------------------------------------
image2mask.py:

Provides a function to obtain mask images based on a given intensity image.  

def image2mask(image_file, roi_image_file='', reverse_pixels_in_roi=False, 
               thresh='otsu', intensity_cutoff=False, cutoff_min=0, 
               cutoff_max=255, cluster_min_size=0, cluster_max_size=0, 
               save_images=False, save_table=False, save_dir='')
               
Input parameters:
image_file - the input image file name
roi_image_file - the roi mask image file name
reverse_pixels_in_roi - True if the interior of the ROI is 0/False
thresh - type of threshold to use, can be a single value, e.g. 'otsu' or a tuple, e.g. ('otsu','kmeans','otsu'), specifying a thresh to apply for each channel (if single value, same thresh will be used for all channels).  The following values are accepted: 'otsu', 'kmeans'.  If given a blank value or any other value, the thresh is set to 0 and all nonzero pixels will be included.
intensity_cutoff - True/False, whether to apply an intensity cutoff to each channel prior to thresholding
cutoff_min - the min intensity value, pixels below this value will be ignored (all channels)
cutoff_max - the max intensity value, pixels above this value will be ignored (all channels)
cluster_min_size - min size of clusters after thresholding, objects below this minimum will be removed (this input can be a single value, e.g. 20 or a tuple, e.g. (20,0,0) where a min size is specified for each color channel separately
cluster_max_size - max size of clusters after thresholding, objects above this maximum will be removed (can be a single value or a tuple).  If set to 0, no max size will be applied.
save_images - If True, a mask is saved for each channel after thresholding (c1_mask.tif, c2_mask.tif, c3_mask.tif) and also a combined mask image (RGB): c123_mask.tif.
save_table - If True, tables are saved (c1_table.csv, c2_table.csv, c3_table.csv) for each channel listing the following properties for each object in the image (after thresholding): 
'area', 'perimeter', 'centroid', 'weighted_centroid', 'major_axis_length', 'minor_axis_length', 'eccentricity', 'orientation', 'min_intensity', 'mean_intensity', 'max_intensity'
save_dir - The directory to save the images/table.

Return value:
A tuple: (mask_image, roi_selection, intensity_image, image_measurements), where 
mask_image = the RGB mask image, each channel is of type bool and has been thresholded, etc. according the input parameters.
roi_selection = the roi mask image.
intensity_image = the original intensity image read from the file
image_measurements = a dict containing summary information for each of the channels of the image.  The values of the dict are tuples, as detailed below.
image_measurements['cluster_total_area'] - (c1-value, c2-value, c3-value)
image_measurements['cluster_count'] - (c1-value, c2-value, c3-value)
image_measurements['cluster_ave_area'] - (c1-value, c2-value, c3-value)
image_measurements['overlap_total_area'] - ((c1+c1-value,c1+c2-value,c1+c3-value), (c2+c1-value,c2+c2-value,c2+c3-value), (c3+c1-value,c3+c2-value,c3+c3-value))
image_measurements['overlap_count'] - ((c1+c1-value,c1+c2-value,c1+c3-value), (c2+c1-value,c2+c2-value,c2+c3-value), (c3+c1-value,c3+c2-value,c3+c3-value))

----------------------------------------------------------------------------------------------------------------------------
Additional scripts included as examples of how to use the module:

run_calculate_IF.py - a script illustrating/testing the IF calculation function.
run_simulate_IF.py - a script illustrating/testing the IF simulation function.
run_simulate_IF_ellipses.py - a script illustrating/testing the IF simulations using ellipses instead of the original clusters of the image.

The input image/ROI used in the example scripts is included in this repository.



