#! /usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:                                                                       #
#                                                                             #
# PURPOSE:  A script to use a KDTree for a computationally efficient          #
#           nearest-neighbour search.                                         #
#           Original created by Jamie Farnes - USyd - 20th September 2012.    #
#                                                                             #
# MODIFIED: 14-May-2015 by C. Purcell                                         #
#                                                                             #
#=============================================================================#
import numpy as np
import time
import scipy
from scipy import spatial

print "Begin using KDTree..."
def do_kdtree(combined_x_y_arrays, points, num_nearest_neighbours):
	
	mytree = scipy.spatial.cKDTree(combined_x_y_arrays)
	# Return the "num_nearest_neighbours" nearest neighbours
	dist, indexes = mytree.query(points, num_nearest_neighbours)
	return dist, indexes

# arrays of the x or y values i.e. [ra1 ra2 ra3 ra4]
#x_array = ra_1; y_array = dec_1

# array of the x/y values in the form i.e. [(ra1,dec1) (ra2,dec2) (ra3,dec3)]
#points = DATA2_coords

x_array = np.array([ 0.6, 30.2, 100.1 ])
y_array = np.array([ 3.7, 15.3, 109.2 ])
points = np.array([(0,0),(20,20),(100.1,109.3)])

# Prepare the data for entry into the KDTree routine
combined_x_y_arrays = np.dstack([x_array.ravel(), y_array.ravel()])[0]
start = time.time()
distances, results = do_kdtree(combined_x_y_arrays, points,1)
end = time.time()
print 'Completed in: ', end-start, "seconds."

# Data comes out of this tree in form (x, y) such that:
# points[0] = Some source coordinates
# x_array[results[0]], y_array[results[0]] = the nearest source to points[0].
# distances[0] = the distance to the nearest source.
