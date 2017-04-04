# -*- coding: utf-8 -*-
"""
#    Copyright (C) 2017  Sarah Keegan
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

import interaction_factor as IF
from image2mask import image2mask
test_colors = False

dir_ = '/Users/sarahkeegan/images/'
image_file_name = dir_ + 'cell-1_1_R_G.tif'
roi_file_name = dir_ + 'cell-1_1_ROI.tif'

ret_images = image2mask(image_file_name, roi_file_name, save_images=True, 
                        save_table=True, save_dir=dir_)
                        
#image2mask image measurements:
print
print "image2mask image measurements:"
for i in range(3):
    print "Channel " + str(i+1) + ":"
    for key in ret_images[3].keys():
        print key + ": " + str(ret_images[3][key][i])
    print

#initialize class
my_IF = IF.interaction_factor(ret_images[0], ret_images[1], ret_images[2])

#create a simulation using the segmented clusters of the image   
simulated_IF = 0.75
my_IF.save_dir=dir_     

#test colors  
if(test_colors):
    my_IF.save_colors = ['red','blue']
    my_IF.simulate_IF(simulated_IF)
    
    my_IF.save_colors = ['magenta','green']
    my_IF.simulate_IF(simulated_IF)
    
    #back to default, run final simulation
    my_IF.save_colors = ['red','green']

my_IF.simulate_IF(simulated_IF) 
print      
print "Finished simulation at R-G IF=" + str(simulated_IF)

#output simulated image measurements:
print
print "Simulated image measurements:"
print "num overlaps, ch" + str(my_IF.nonref_color+1) + ", ch" + str(my_IF.ref_color+1) + ": " + str(my_IF.sim_num_ov_clusters)
print "total area overlap clusters: " + str(my_IF.sim_area_ov_clusters)
print "% overlapping ch" + str(my_IF.ref_color+1) + " clusters: " + str(my_IF.sim_num_ref_ov_clusters/float(my_IF.num_clusters[my_IF.ref_color]))
    
#calculate IF of the simluated image just created
sim_image_file_name = dir_ + my_IF.save_colors[0] + '-' + my_IF.save_colors[1] + '_IF_' + str(simulated_IF) + '_sim_mask.tif'
ret_images = image2mask(sim_image_file_name, roi_file_name)

#image2mask image measurements:
print
print "image2mask image measurements (of image simulated at R-G IF="+str(simulated_IF)+"):"
for i in range(3):
    print "Channel " + str(i+1) + ":"
    for key in ret_images[3].keys():
        print key + ": " + str(ret_images[3][key][i])
    print
    
my_IF2 = IF.interaction_factor(ret_images[0], ret_images[1], ret_images[2])
ret_vals = my_IF2.calculate_IF()
print
print "Calculated R-G IF of the simulated image (R-G IF "+str(simulated_IF)+") is: " + str(ret_vals[0])
    
#output 'original' image measurements, as a check:
print
print "Original image measurements (of image simulated at R-G IF="+str(simulated_IF)+"):"
for i in range(3):
    print "total area of clusters, ch" + str(i+1) + ": " + str(my_IF2.orig_area_clusters[i])
print "num overlaps, ch" + str(my_IF2.nonref_color+1) + ", ch" + str(my_IF2.ref_color+1) + ": " + str(my_IF2.orig_num_ov_clusters)
print "total area overlap clusters: " + str(my_IF2.orig_area_ov_clusters)
print "% overlapping ch" + str(my_IF2.ref_color+1) + " clusters: " + str(my_IF2.orig_num_ref_ov_clusters/float(my_IF2.num_clusters[my_IF2.ref_color]))
    
#calculate reverse IF for the image
my_IF2.nonref_color = 1
my_IF2.ref_color = 0
ret_vals = my_IF2.calculate_IF()
print
print "Calculated G-R IF of the simulated image (R-G IF "+str(simulated_IF)+") is: " + str(ret_vals[0])

#output 'original' image measurements (again, only difference should be % overlapping since we changed the ref channel):
print
print "Original image measurements (of image simulated at R-G IF="+str(simulated_IF)+"):"
for i in range(3):
    print "total area of clusters, ch" + str(i+1) + ": " + str(my_IF2.orig_area_clusters[i])
print "num overlaps, ch" + str(my_IF2.nonref_color+1) + ", ch" + str(my_IF2.ref_color+1) + ": " + str(my_IF2.orig_num_ov_clusters)
print "total area overlap clusters: " + str(my_IF2.orig_area_ov_clusters)
print "% overlapping ch" + str(my_IF2.ref_color+1) + " clusters: " + str(my_IF2.orig_num_ref_ov_clusters/float(my_IF2.num_clusters[my_IF2.ref_color]))
    

     