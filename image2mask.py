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

from scipy import ndimage
import numpy as np
from skimage import io, exposure, measure
from skimage.filters import threshold_otsu
from sklearn.cluster import k_means
import warnings
import pandas as pd

warnings.filterwarnings('ignore', message='.+ is a low contrast image', category=UserWarning)

labeling_structure = [[1,1,1],[1,1,1],[1,1,1]]
#[[0,1,0],[1,1,1],[0,1,0]]

image_prop_list = ['area', 'perimeter', 'centroid', 'weighted_centroid', 
                   'major_axis_length', 'minor_axis_length', 'eccentricity', 
                   'orientation', 'min_intensity', 'mean_intensity', 'max_intensity']

#Input: RGB image file, ROI file
#Output: C1 mask, C2 mask, C3 mask, ROI selection OR bbox limits
#Parameters = threshold method

def image2mask(image_file, roi_image_file='', reverse_pixels_in_roi=False, 
               thresh='otsu', intensity_cutoff=False, cutoff_min=0, 
               cutoff_max=255, cluster_min_size=0, cluster_max_size=0, 
               save_images=False, save_table=False, save_dir=''):
    
    if(type(thresh) != type([])):
        thresh = [thresh,thresh,thresh]
        
    if(type(cluster_min_size) != type([])):
        cluster_min_size = [cluster_min_size,cluster_min_size,cluster_min_size]
        
    if(type(cluster_max_size) != type([])):
        cluster_max_size = [cluster_max_size,cluster_max_size,cluster_max_size]
    
    #roi to mask
    if(roi_image_file != ''):
        roi_image = io.imread(roi_image_file)
        if len(roi_image.shape) == 3:
            roi_image = roi_image[:,:,1] 
            
        if(reverse_pixels_in_roi):
            roi_mask = roi_image > 0            
            #roi_mask is TRUE outside of ROI
        else:
            roi_mask = roi_image == 0            
            #roi_mask is TRUE outside of ROI
        roi_selection = ~np.array(roi_mask) 
        #roi_selection is TRUE inside of ROI
    else:
        roi_mask = []
        roi_selection = []
    
    #split to 3 colors
    rgb_image = io.imread(image_file)
    split_images = []

    if(roi_image_file == ''):
        roi_mask = np.zeros_like(rgb_image[:,:,0])
        roi_mask = roi_mask.astype(dtype='bool')
        roi_selection = ~np.array(roi_mask)

    for color_i in range(3):
        
        cur_image = rgb_image[:,:,color_i]
        
        cur_image[roi_mask] = 0 #set all pixels outside roi to 0
            
        #intensity cutoff
        if(intensity_cutoff):
            min_in_range = cutoff_min
            max_in_range = cutoff_max
            cur_image = exposure.rescale_intensity(cur_image, 
                                                   in_range=(min_in_range,max_in_range), 
                                                   out_range=(0,255))
        #threshold
        if(cur_image.max() == 0): cur_th = 0
        else:
            if(type(thresh[color_i]) == type(0) or type(thresh[color_i]) == type(0.0)):
                cur_th = thresh[color_i] #manual threshold option
            else:
                if thresh[color_i] == 'otsu':
                    cur_th = threshold_otsu(cur_image[roi_selection == True])
                elif thresh[color_i] == 'kmeans':
                    pixels = np.expand_dims(cur_image[roi_selection == True].flatten(),axis=1)
                    ret = k_means(pixels, 3)
                    intensity_cutoffs = sorted(ret[0])
                    cur_th = intensity_cutoffs[2][0]
                else:
                    cur_th = 0
        
        cur_mask = cur_image > cur_th
            
        if(cluster_min_size[color_i] > 0 or cluster_max_size[color_i] > 0):
            #remove clusters outside of constraints from the mask
            #find measurements for image, area and count of clusters
            labeled_clusters, num_clusters = ndimage.label(cur_mask, structure=labeling_structure)
            props = measure.regionprops(labeled_clusters)
            for region in props:
                if(region['area'] < cluster_min_size[color_i] or (cluster_max_size[color_i] > 0 and region['area'] > cluster_max_size)):
                    #remove region
                    coords = region['coords']
                    for cur_coord in coords:
                        cur_mask[cur_coord[0],cur_coord[1]] = False
                    
        #save mask and image
        if(save_images):
            io.imsave(save_dir + '/c' + str(color_i+1) + '_mask.tif', 255*cur_mask.astype('uint8')) 
            
        split_images.append([cur_image,cur_mask])    
        
    #find measurements for image, area and count of clusters
    image_measurements = {}
    for color in range(3):
        if(color==0):
            image_measurements['cluster_total_area'] = []
            image_measurements['cluster_count'] = []
            image_measurements['cluster_ave_area'] = []
            image_measurements['overlap_total_area'] = []
            image_measurements['overlap_count'] = []
        image_measurements['cluster_total_area'].append(split_images[color][1].sum())
        
        labeled_clusters, num_clusters = ndimage.label(split_images[color][1], structure=labeling_structure)
        image_measurements['cluster_count'].append(num_clusters)
        props = measure.regionprops(labeled_clusters, split_images[color][0])
        ave_area = 0
        for region in props:
            ave_area += region.area 
        if(len(props) > 0): ave_area /= float(len(props))
        image_measurements['cluster_ave_area'].append(ave_area)
        
        ov_area_arr = []
        ov_count_arr = []
        for ov_color in range(3):
            if(ov_color != color):
                and_area = np.logical_and(split_images[color][1],split_images[ov_color][1])
                ov_area_arr.append(and_area.sum())
                
                labeled_clusters, num_clusters = ndimage.label(and_area, structure=labeling_structure)
                ov_count_arr.append(num_clusters)
        
        image_measurements['overlap_total_area'].append(ov_area_arr)
        image_measurements['overlap_count'].append(ov_count_arr)
        
        if(save_table):
            #save entire particles table for each set of clusters
            prop_table = {}
            for prop in image_prop_list:
                prop_table[prop] = []
                for region in props:
                    prop_table[prop].append(region[prop])
                
            prop_df = pd.DataFrame(prop_table)
            prop_df.to_csv(save_dir + '/c' + str(color+1) + '_table.csv')
    
    intensity_image = np.dstack((split_images[0][0],split_images[1][0], split_images[2][0]))
    mask_image = np.dstack((split_images[0][1], split_images[1][1], split_images[2][1]))
    
    #save rgb mask and intensity image
    if(save_images):
        io.imsave(save_dir + '/c123_mask.tif', mask_image.astype('uint8')*255)
        
    return (mask_image, roi_selection, intensity_image, image_measurements)



