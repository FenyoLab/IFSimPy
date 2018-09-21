# -*- coding: utf-8 -*-
"""
#    Copyright (C) 2017  Sarah Keegan and Keria Bermudez-Hernandez
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from skimage import io, measure
from skimage.draw import ellipse_perimeter
from skimage.transform import rotate
from scipy import ndimage
import random

"""
This class will calculate the "Interaction Factor (IF)" for an image.

"""
class interaction_factor:
    
    """
    Initialize class: a class instance should be created for an RGB mask (segmented) image
    The following class variables are set or initialized to default values:
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
    
    """
    def __init__(self, rgb_mask_image, roi_mask=[], intensity_image=[]):
        
        # Image parameters
        self.labeling_structure = [[1,1,1],[1,1,1],[1,1,1]]   
        self.cluster_masks = []
        self.num_clusters = []
        self.cluster_props = []
        for i in range(3):
            self.cluster_masks.append(rgb_mask_image[:,:,i])
            labeled_clusters, num_clusters = ndimage.label(self.cluster_masks[i], 
                                                           structure=self.labeling_structure)
            self.num_clusters.append(num_clusters)
            self.cluster_props.append(measure.regionprops(labeled_clusters))
            
        #set up bounding box of ROI
        self.roi_mask = roi_mask
        if(roi_mask == []):
            self.use_roi = False
            self.area_roi = len(rgb_mask_image) * len(rgb_mask_image[0]) #area ROI is entire image area
            self.roi_min_row = 0
            self.roi_min_col = 0
            self.roi_max_row = len(rgb_mask_image)-1
            self.roi_max_col = len(rgb_mask_image[0])-1
        else:
            self.use_roi = True
            labeled_clusters, num_clusters = ndimage.label(self.roi_mask, 
                                                           structure=self.labeling_structure)
            props = measure.regionprops(labeled_clusters)
            self.area_roi = props[0].area
            bounding_box = props[0].bbox
            self.roi_min_row = int(bounding_box[0])
            self.roi_min_col = int(bounding_box[1])
            self.roi_max_row = int(bounding_box[2])
            self.roi_max_col = int(bounding_box[3])
            
        self.intensity_image = intensity_image

        # IF Parameters
        self.nonref_color=0
        self.ref_color=1
        self.move_nonref=True
        
        # Parameters used in calculate_IF
        self.plot_IF_curve=False
        self.save_random_images=False
        self.num_random_sims=50
        
        # Parameters used in simulate_IF
        self.save_image=True
        self.save_colors=['red','green'] #[non-ref, ref], choices: red, green, magenta (should only be paired with green)
        self.image_file_prefix=''
        self.use_ellipses=False
        self.nonref_size=1
        self.ref_size=1
        self.nonref_number=None
        self.ref_number=None
        self.rotate_ref=True
        
        # Directory for saving output of either function
        self.save_dir=''

        # Setting to allow nearby clusters to be consdiered as 'overlap' for IF calculation
        # expand_pixels indicates how much of pixel border around cluster to consider as 'overlapping'
        self.allow_nearby_clusters = False
        self.expand_pixels = 1
                 
    """ 
    Calculates measurements for the input (orig) image 
    Rerun this if changing ref/nonref clusters
    Calculates the following measurements:
    --------------------------------------
    orig_area_clusters      : array, sum of total area of clusters, each channel
    orig_overlap_mask       : overlap mask for ref+nonref colors 
    orig_num_ov_clusters    : number of overlaps in overlap mask
    orig_area_ov_clusters   : sum of total area of clusters in overlap mask
    orig_num_ref_ov_clusters: number of reference clusters that overlap non-ref clusters
    """
    def calc_orig_image_measurements(self):
    
        self.orig_area_clusters = []
        for i in range(3):
            self.orig_area_clusters.append(self.cluster_masks[i].sum())
        
        #provides overlap information for ref/nonref channel clusters
        self.orig_overlap_mask = np.logical_and(self.cluster_masks[self.nonref_color], self.cluster_masks[self.ref_color])
        labeled_clusters, self.orig_num_ov_clusters = ndimage.label(self.orig_overlap_mask, structure=self.labeling_structure)
        self.orig_area_ov_clusters = self.orig_overlap_mask.sum()
        
        #get count reference clusters that overlap
        labeled_clusters, num_clusters = ndimage.label(self.cluster_masks[self.ref_color], structure=self.labeling_structure)
        props = measure.regionprops(labeled_clusters)
        self.orig_num_ref_ov_clusters = 0
        height = len(self.cluster_masks[self.nonref_color])
        width = len(self.cluster_masks[self.nonref_color][0]) #check!
        for region in props:
            coords = region.coords
            overlap_nonref = False
            for i in np.arange(0,len(coords)):
                y = coords[i][0]
                x = coords[i][1]
                if (self.allow_nearby_clusters):
                    # first make sure it's a border pixel
                    border_pixel = False
                    for x_add_px in [-1, 0, 1]:
                        for y_add_px in [-1, 0, 1]:
                            if (y + y_add_px >= height or x + x_add_px >= width or y + y_add_px < 0 or x + x_add_px < 0 or
                                    (not self.cluster_masks[self.ref_color][y + y_add_px, x + x_add_px])):
                                border_pixel = True
                    # then look for overlap
                    if (border_pixel):
                        for x_add_px in range(-1 * self.expand_pixels, self.expand_pixels + 1, 1):
                            for y_add_px in range(-1 * self.expand_pixels, self.expand_pixels + 1, 1):
                                if (y + y_add_px < height and x + x_add_px < width and y + y_add_px >= 0 and x + x_add_px >= 0 and
                                        self.cluster_masks[self.nonref_color][y + y_add_px, x + x_add_px]):
                                    overlap_nonref = True
                                    break
                    elif(self.cluster_masks[self.nonref_color][y, x]): #(if not on border, only look at the single position)
                        overlap_nonref=True
                        break
                elif(self.cluster_masks[self.nonref_color][y, x]):
                    overlap_nonref = True
                    break

            if(overlap_nonref): 
                self.orig_num_ref_ov_clusters += 1

    def draw_outline_nearby_clusters(self):
        #draws expanded outlines around nonref clusters in rgb image mask,
        #based on expand_pixels setting, for use as a display image when using allow_nearby_clusters setting
        #returns rgb image with outlines and also saves to save_dir

        outline_image = np.zeros((len(self.cluster_masks[self.nonref_color]), len(self.cluster_masks[self.nonref_color][0]), 3), dtype='uint8')
        self._ch2color(outline_image, self.cluster_masks[self.nonref_color], self.save_colors[0])
        self._ch2color(outline_image, self.cluster_masks[self.ref_color], self.save_colors[1])
        outline_image = 255 * outline_image.astype('uint8')

        labeled_clusters, num_clusters = ndimage.label(self.cluster_masks[self.ref_color], structure=self.labeling_structure)
        props = measure.regionprops(labeled_clusters)
        height = len(self.cluster_masks[self.nonref_color])
        width = len(self.cluster_masks[self.nonref_color][0])

        for region in props:
            coords = region.coords
            for i in np.arange(0, len(coords)):
                y = coords[i][0]
                x = coords[i][1]

                # first make sure it's a border pixel
                border_pixel=False
                for x_add_px in [-1,0,1]:
                    for y_add_px in [-1,0,1]:
                        if (y + y_add_px >= height or x + x_add_px >= width or y + y_add_px < 0 or x + x_add_px < 0 or
                                (not self.cluster_masks[self.ref_color][y + y_add_px, x + x_add_px])):
                            border_pixel=True

                #then draw outline
                if(border_pixel):
                    for x_add_px in range(-1*self.expand_pixels,self.expand_pixels + 1, 1):
                        for y_add_px in range(-1*self.expand_pixels,self.expand_pixels + 1, 1):
                            if (y + y_add_px < height and x + x_add_px < width and y + y_add_px >=0 and x + x_add_px >= 0 and
                                    (not self.cluster_masks[self.ref_color][y + y_add_px, x + x_add_px]) and
                                    (not self.cluster_masks[self.nonref_color][y + y_add_px, x + x_add_px])):
                                outline_image[y + y_add_px, x + x_add_px] = [255, 255, 255]

        #save
        if(self.save_dir != ''):
            io.imsave(self.save_dir + '/expanded_clusters-orig_image.tif', outline_image)

        return outline_image

    """ 
    Calculates the IF for an image 
    - returns tuple containing: (IF, p-value, % of ref clusters overlapping nonref clusters for original image)
    - plots the IF equation and saves to the save_dir if plot_IF_curve is set to True
    """
    def calculate_IF(self):
        
        #run the random simulations
        x_data = [] # percentage of overlap
        cluster_ov_stats = [] #each cluster, whether it overlaps 
        
        save_image_value = self.save_image #set this back at the end
        file_prefix_value = self.image_file_prefix
        self.save_image = self.save_random_images
        
        for sim_i in range(self.num_random_sims):
            if(self.save_random_images): self.image_file_prefix = str(sim_i+1)
            self.simulate_IF() #run sim for IF==0

            print(str(sim_i) + ' ' + str(self.sim_num_ref_ov_clusters) + ' ' + str(float(self.sim_num_ref_ov_clusters) / self.num_clusters[self.ref_color]))

            x_data.append(float(self.sim_num_ref_ov_clusters) /
                                self.num_clusters[self.ref_color])
            if(sim_i == 0):
                orig_meas = (float(self.orig_num_ref_ov_clusters) / 
                                   self.num_clusters[self.ref_color])
                               
            cluster_ov_stats.append(self.ref_cluster_ov)
        
        #print "FINISHED SIM IMAGES"
        
        #Use IF equation to calculate IF
        #for each cluster, k, get probability of each cluster to overlap, i.e. P0k
        cluster_list = cluster_ov_stats[0].keys()
        ov_prob = []
        for k in cluster_list:
            ov_count = 0
            for sim_i in cluster_ov_stats:
                ov_count += sim_i[k]
            ov_prob.append(ov_count/float(len(cluster_ov_stats)))
            
        #for each IF, 0 to 1, plug in to formula and see how close it approximates
        #the experimental % of overlap ('P')
        calcIF = 0
        IF_increment = 0.01
        min_calcIF = 0
        diff = 0
        x = np.arange(min_calcIF,1+IF_increment,IF_increment)
        y = []
        for IF in x:
            cur_m = self._IF_equation(IF, ov_prob)
            y.append(cur_m)
            cur_diff = abs(cur_m - orig_meas)
            if(IF == min_calcIF or cur_diff < diff):
                calcIF = IF
                diff = cur_diff
            if(not self.plot_IF_curve or self.save_dir == ''):
                #can stop, found best IF, but if we want to plot the curve, keep going
                if(cur_diff > diff):
                    break
        count_p = 0
        for val in x_data:
            if(val > orig_meas):
                count_p += 1  
        p_val = count_p / float(len(x_data)) 
             
        if(self.plot_IF_curve and self.save_dir != ''):
            fig1 = plt.figure()
            ax = fig1.add_subplot(111)
            ax.axhline(y=orig_meas)
            ax.axvline(x=calcIF, color='red')
            ax.set_title('IF = ' + str(calcIF) + ', %-overlap = ' + str(100*round(orig_meas,3)))
            ax.plot(x, y, '-')
            ax.set_xlim(0,1)
            ax.set_ylim(0,1)
            fig1.savefig(self.save_dir + '/IF_Equation')
            ax.cla()
            fig1.clf()
            
        self.save_image = save_image_value
        self.image_file_prefix = file_prefix_value
        
        return (calcIF, p_val, orig_meas)
    
    """
    Performs a simulation at given IF level
    Simulated image will be output to a file if self.save_image==True and self.save_dir!=''
    Calculates/creates the following measurements/images for the simulated image:
    --------------------------------------------------------------
    sim_num_ref_ov_clusters: number of reference overlapping clusters in the simulation
    sim_overlap_mask       : overlap mask for ref+nonref clusters in the simulation
    sim_num_ov_clusters    : number of overlaps in overlap mask (simulation)
    sim_area_ov_clusters   : sum of total area of clusters in overlaop mask (simulation)
    simulated_images       : array, simulated image for each channel, includes intensities, 
                           : if intensity image is not empty (otherwise values 
                             are either 0 or 255)
    simulated_masks        : array, simulated image for each channel,  boolean values
    ref_cluster_ov         : dict, T/F whether ref cluster with given label overlaps
                           : a nonref cluster (label is from cluster_props[ref_color])
    """
    def simulate_IF(self, IF=0):
        
        #NOTE: if numbers/ellipses are changed here, can't be changed back in a subsequent run
        
        #if given a new number of clusters (rather than all clusters in image),
        #randomly select that number of clusters from the prop list and use these for the simulation
        #(this option is used when making simulations with different numbers of clusters to test
        #that IF is not depandant on density)      
        if(self.nonref_number != None):
            self._adjust_cluster_numbers(self.nonref_number, self.nonref_color)
        if(self.ref_number != None):
            self._adjust_cluster_numbers(self.ref_number, self.ref_color)
            
        #use ellipses instead of clusters for the simulation, size is specified
        #(this option is to change from clusters to ellipses in order to make simulations
        #with different size ellipses to test that IF is not dependant on density)
        if self.use_ellipses:
            self._ellipses(self.nonref_color, self.nonref_size)
            self._ellipses(self.ref_color, self.ref_size)
            
        #can't produce precision higher than 0.01
        IF = round(IF, 2)
        
        #beginning simulation
        self.simulated_images = [[],[],[]]
        
        #(1) placement of non-reference color
        if(self.move_nonref):

            #generating random image, restricted to bbox of the ROI
            self.simulated_images[self.nonref_color] = np.zeros(self.cluster_masks[self.nonref_color].shape,dtype='uint8')
            for region in self.cluster_props[self.nonref_color]:
                bounding_box = region.bbox
                coords = region.coords
                min_row = int(bounding_box[0])
                min_col = int(bounding_box[1])
                max_row = int(bounding_box[2])
                max_col = int(bounding_box[3])
                region_height = max_row-min_row
                region_width = max_col-min_col
                cluster_at_origin = coords - np.array([min_row,min_col])

                max_num_tries = 10000
                num_tries = 0
                while num_tries < max_num_tries:
                    #generating random number and new coordinates 
                    random_row = random.randint(self.roi_min_row, self.roi_max_row-region_height)
                    random_col = random.randint(self.roi_min_col, self.roi_max_col-region_width)
                    new_coords = cluster_at_origin + np.array([random_row,random_col])
                    
                    placement_failed = False
                    #if placement_failed, try again
                    if ((self.use_roi and self._outside_roi(new_coords)) or self._overlapping(new_coords, self.nonref_color)):
                       placement_failed = True
                       num_tries += 1
                    else:
                        break
                        
                if(placement_failed):
                    print("Warning: placement of non-reference cluster failed: max_num_tries exceeded.")
                    #print warning to error log since exceeded max_tries
                else:
                    #place cluster
                    for i in np.arange(0,len(new_coords)):
                        row_coords = new_coords[i][0]
                        column_coords = new_coords[i][1]
                        orig_row_coords = coords[i][0]
                        orig_column_coords = coords[i][1]
                        if(not self.use_ellipses and len(self.intensity_image) > 0):
                            self.simulated_images[self.nonref_color][row_coords,column_coords] = (
                                self.intensity_image[orig_row_coords,orig_column_coords][self.nonref_color])
                        else:
                            self.simulated_images[self.nonref_color][row_coords,column_coords] = 255
                            
        else:
            if(self.use_ellipses):
                self.simulated_images[self.nonref_color] = np.zeros(self.cluster_masks[self.nonref_color].shape,
                                                                dtype='uint8')
                for region in self.cluster_props[self.nonref_color]: #these are the ellipses
                    coords = region.coords
                    #self.simulated_images[self.nonref_color][coords[i][0],coords[i][1]] = 255
                    for i in len(coords):
                        self.simulated_images[self.nonref_color][coords[i][0], coords[i][1]] = 255
            else:
                #create 'simulated image' from mask image, which is just the mask image changed to intensity image
                self.simulated_images[self.nonref_color] = self.intensity_image[:,:,self.nonref_color] * self.cluster_masks[self.nonref_color].astype('uint8')
                 
        #(2) placement of reference color
        self.simulated_images[self.ref_color] = np.zeros(self.cluster_masks[self.ref_color].shape,
                                                         dtype='uint8')
                                                          
        #the number of ref clusters that overlap a non-ref cluster
        self.sim_num_ref_ov_clusters = 0
        
        #for each ref cluster, stores whether it is overlapping a nonref cluster (True/False)
        self.ref_cluster_ov = {}

        height = len(self.simulated_images[self.nonref_color])
        width = len(self.simulated_images[self.nonref_color][0])  # check!
        for region in self.cluster_props[self.ref_color]:
            label = region.label
            bounding_box = region.bbox
            coords = region.coords
            maj_axis = region.major_axis_length
            if(maj_axis == 0):
                maj_axis = 1
            min_row = int(bounding_box[0])
            min_col = int(bounding_box[1])
            max_row = int(bounding_box[2])
            max_col = int(bounding_box[3])
            region_height = max_row-min_row
            region_width = max_col-min_col
            cluster_at_origin = coords - np.array([min_row,min_col])

            if(self.rotate_ref):
                # cluster rotation- (1) choose random angle, (2) transform
                if(len(region.image) >=2 and len(region.image[0]) >= 2):
                    angle = random.randint(0, 360)
                    cluster_rot_img = rotate(region.image.astype('uint8'), angle, resize=True)
                    cluster_rot_img = cluster_rot_img>0

                    #rotation leaves border, get rid of...
                    l, n = ndimage.label(cluster_rot_img, structure=self.labeling_structure)
                    p = measure.regionprops(l)

                    l, n = ndimage.label(p[0].image, structure=self.labeling_structure)
                    p = measure.regionprops(l)
                    cluster_rot_coords = p[0].coords

                    region_height = int(p[0].bbox[2]) - int(p[0].bbox[0])
                    region_width = int(p[0].bbox[3]) - int(p[0].bbox[1])

                    #io.imsave('/Users/sarahkeegan/fenyolab/temp/orig.tif', region.image.astype('uint8')*255)
                    #io.imsave('/Users/sarahkeegan/fenyolab/temp/rot.tif', p[0].image.astype('uint8')*255)

                else:
                    cluster_rot_coords = cluster_at_origin
            
            max_num_tries = 10000
            num_tries = 0
            while num_tries < max_num_tries:
                #generating random number and new coordinates 
                random_row = random.randint(self.roi_min_row, self.roi_max_row-region_height)
                random_col = random.randint(self.roi_min_col, self.roi_max_col-region_width)
                if(self.rotate_ref):
                    new_coords = cluster_rot_coords + np.array([random_row,random_col])
                else:
                    new_coords = cluster_at_origin + np.array([random_row, random_col])
                
                placement_failed = False
                if ((self.use_roi and self._outside_roi(new_coords)) or self._overlapping(new_coords, self.ref_color)):
                       placement_failed = True
                       num_tries += 1 #re-try placement
                else:
                    #place cluster based on IF
                    #uniform (IF == 0) or with increased probability of overlap (IF > 0)
                    overlap_nonref = False
                    for i in np.arange(0,len(new_coords)):
                        y = new_coords[i][0]
                        x = new_coords[i][1]
                        if (self.allow_nearby_clusters):
                            # first check if it's a border pixel
                            border_pixel = False
                            for x_add_px in [-1, 0, 1]:
                                for y_add_px in [-1, 0, 1]:
                                    if (y + y_add_px >= height or x + x_add_px >= width or y + y_add_px < 0 or x + x_add_px < 0 or
                                            (self.simulated_images[self.ref_color][y + y_add_px, x + x_add_px] == 0)):
                                        border_pixel = True
                            # then look for overlap
                            if (border_pixel):
                                for x_add_px in range(-1 * self.expand_pixels, self.expand_pixels + 1, 1):
                                    for y_add_px in range(-1 * self.expand_pixels, self.expand_pixels + 1, 1):
                                        if (y + y_add_px < height and x + x_add_px < width and y + y_add_px >= 0 and x + x_add_px >= 0 and
                                                (self.simulated_images[self.nonref_color][y + y_add_px, x + x_add_px] != 0)):
                                            overlap_nonref = True
                                            break
                            elif(self.simulated_images[self.nonref_color][y,x] != 0):
                                overlap_nonref = True
                                break
                        elif(self.simulated_images[self.nonref_color][y,x] != 0):
                            overlap_nonref = True
                            break
                    if(IF != 0):
                        if(not overlap_nonref):
                            #no overlap with non-ref cluster, decide whether to re-try cluster placement based on IF,
                            #flip "bias coin": e.g. if IF == 0.75, probability of re-try is .75 (and keep is .25)
                            random_number = round(random.random(),2) #we do not expect to predict IF to more precision than .01
                            if(random_number <= IF): 
                                num_tries += 1 #re-try placement
                            else:
                                break
                        else: 
                            break
                    else:
                        break
            if(placement_failed):
                print("Warning: placement of reference cluster failed: max_num_tries exceeded.")
                #print warning to error log since exceeded max_tries
            else:
                #place cluster
                for i in np.arange(0,len(new_coords)):
                    row_coords = new_coords[i][0]
                    column_coords = new_coords[i][1]
                    if(not self.rotate_ref):
                        orig_row_coords = coords[i][0]
                        orig_column_coords = coords[i][1]
                    if(not self.use_ellipses and len(self.intensity_image) > 0 and not self.rotate_ref):
                        self.simulated_images[self.ref_color][row_coords,column_coords] = (
                            self.intensity_image[orig_row_coords,orig_column_coords][self.ref_color])
                    else:
                        self.simulated_images[self.ref_color][row_coords,column_coords] = 255
                if(overlap_nonref): 
                    self.sim_num_ref_ov_clusters += 1
                    self.ref_cluster_ov[label] = True
                else:
                    self.ref_cluster_ov[label] = False
                    
        self.simulated_masks = [[],[],[]]
        self.simulated_masks[self.nonref_color] = self.simulated_images[self.nonref_color] > 0
        self.simulated_masks[self.ref_color] = self.simulated_images[self.ref_color] > 0
                    
        #(3) save results
                    
        #calculate sim image measurements
        self.sim_overlap_mask = np.logical_and(self.simulated_masks[self.nonref_color], self.simulated_masks[self.ref_color])
        labeled_clusters, self.sim_num_ov_clusters = ndimage.label(self.sim_overlap_mask, structure=self.labeling_structure)
        self.sim_area_ov_clusters = self.sim_overlap_mask.sum()
        
        #calcuate orig image measurements 
        #(may have already been done but we will redo it here in case ref/nonref was changed before fn call)
        self.calc_orig_image_measurements()
                    
        #save simulated image as RGB
        if(self.save_image):
            if(self.image_file_prefix == ''):
                base_filename = self.save_dir + '/' 
            else:
                base_filename = self.save_dir + '/' + self.image_file_prefix + '_'
            
            #save RGB mask/image
            rgb = np.zeros((len(self.simulated_masks[self.nonref_color]), len(self.simulated_masks[self.nonref_color][0]), 3), dtype='uint8')
            self._ch2color(rgb, self.simulated_masks[self.nonref_color], self.save_colors[0])  
            self._ch2color(rgb, self.simulated_masks[self.ref_color], self.save_colors[1])  
            io.imsave(base_filename + self.save_colors[0] + '-' + self.save_colors[1] + '_IF_' + str(IF) + '_sim_mask.tif', 255*rgb.astype('uint8'))
            
            #save intensity images, if given
            if(not self.use_ellipses and len(self.intensity_image) > 0):
                rgb = np.zeros((len(self.simulated_masks[self.nonref_color]), len(self.simulated_masks[self.nonref_color][0]), 3), dtype='uint8')
                self._ch2color(rgb, self.simulated_images[self.nonref_color], self.save_colors[0])  
                self._ch2color(rgb, self.simulated_images[self.ref_color], self.save_colors[1])                  
                io.imsave(base_filename + self.save_colors[0] + '-' + self.save_colors[1] + '_IF_' + str(IF) + '_sim_image.tif', rgb)
                
    ### PRIVATE FUNCTIONS utilized by calculate_IF ###
    """
    IF equation
    This equation relates the IF to the % of overlapping clusters in an image
    The only other parameter is the probability of overlap of each cluster in 
    the random case (which is determined by random simulations)
    """
    def _IF_equation(self, IF, ov_prob):
        #IF equation: input is IF and the probability of overlap of each cluster at IF=0
        #output is the % of overlapping clusters corresponding to a theoretical 'experimental' image
        if(IF == 1): return 1
        p = 0
        for p0 in ov_prob:
            pk = p0/(1-(1-p0)*IF)
            p += pk
        p = p/len(ov_prob)
        return p
        
    ### PRIVATE FUNCTIONS utilized by simulate_IF ###
        
    def _adjust_cluster_numbers(self, num, color):
        new_cluster_props = []
        for i in range(num):
            random_index = random.randint(0,len(self.cluster_props[color])-1)
            prop = self.cluster_props[color][random_index]
            new_cluster_props.append(prop)
        self.cluster_props[color] = new_cluster_props
        self.num_clusters[color] = num
        
    def _ellipses(self, color, size): 
        new_clusters = []
        image_shape = self.cluster_masks[color].shape
        for cluster in self.cluster_props[color]:
            test_mask = np.zeros(image_shape, bool)
        
            centroid = cluster.centroid
            ycenter = int(round(centroid[0]))
            xcenter = int(round(centroid[1]))
            
            major_axis = int(round(cluster.major_axis_length/2 * size)) 
            minor_axis = int(round(cluster.minor_axis_length/2 * size))
            orientation = cluster.orientation
        
            if (major_axis < 1 or minor_axis < 1):
                major_axis = 1
                minor_axis = 1
               
            rr, cc  = ellipse_perimeter(ycenter,xcenter,major_axis,minor_axis,orientation,image_shape)
            test_mask[rr,cc] = True
            test_mask  = ndimage.morphology.binary_fill_holes(test_mask)
            test_mask = np.where(test_mask == True, 1, 0)
            new_cluster = measure.regionprops(test_mask)
            new_clusters += new_cluster
        self.cluster_props[color] = new_clusters
        
    def _outside_roi(self, coords_cluster):
        outside_roi = False
        for i in np.arange(0,len(coords_cluster)):
            row_coords = coords_cluster[i][0]
            column_coords = coords_cluster[i][1]
            # row_coords > len(self.roi_mask) or column_coords > len(self.roi_mask[0]) or
            if(self.roi_mask[row_coords,column_coords] == False):
                outside_roi = True
                break
        return outside_roi
    
    def _overlapping(self, coords_cluster, color):
        # expanding coords around edge so cluster is not overlapping OR touching
        test_mask = self.simulated_images[color]
        mask_with_cluster = np.zeros_like(test_mask, dtype = np.uint8)
        extra_coords = []
        for coord in coords_cluster:
            allowed_r = []
            allowed_c = []
            if(coord[0] > 0):
                allowed_r.append(coord[0]-1)
            if(coord[0] < len(mask_with_cluster)-1):
                allowed_r.append(coord[0]+1)
            if(coord[1] > 0):
                allowed_c.append(coord[1]-1)
            if(coord[1] < len(mask_with_cluster[1])-1):
                allowed_c.append(coord[1]+1)
            for c in allowed_c:
                extra_coords.append([coord[0],c])
            for r in allowed_r:
                extra_coords.append([r,coord[1]])
            for r in allowed_r:
                for c in allowed_c:
                    extra_coords.append([r,c])
        coords_cluster = np.concatenate((coords_cluster,extra_coords),axis=0)
        rr_cluster  = coords_cluster[:,0]   
        cc_cluster = coords_cluster[:,1]
        mask_with_cluster[rr_cluster,cc_cluster] = 1
        non_zero_vals = test_mask[mask_with_cluster > 0].sum()
        overlapping  = non_zero_vals > 0
        return overlapping
        
    def _ch2color(self, rgb_image, image, color):
        if(color == 'red'):
            rgb_image[:,:,0] = image
        if(color == 'green'):
            rgb_image[:,:,1] = image
        if(color == 'blue'):
            rgb_image[:,:,2] = image
        if(color == 'magenta'):
            rgb_image[:,:,0] = image
            rgb_image[:,:,2] = image
            
    
            
        
        
    
    
    
    
    
        
    