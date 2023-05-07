import numpy as np
from astroquery.ned import Ned
import astropy.units as u 
from astropy.coordinates import SkyCoord

def find_angle_wrt_north(coord, galaxy_disk_coord):
    '''
    Takes in a coordinate (usually the one dimensional center of mass at a given distance out from 
    the disk) along with the coordinate of the galaxy disk centroid and calculates the position angle 
    of the line formed between those two points with respect to north.  
    
    Parameters
    ----------
    coord: int array
        x and y coordinates of the point you wish to use in the line for the angle calculation  
    galaxy_disk_coord: int array
        x and y coordinates of the galaxy disk's centroid 

    Returns
    -------
    tail angle: float
        the measured position angle of the line joining the two inputted points with respect 
        to a line pointing north in units of degrees
    '''

    # Finding the difference in x and y positions between the two points to find the absolute value of the angle
    x_diff = coord[0] - galaxy_disk_coord[0]
    y_diff = coord[1] - galaxy_disk_coord[1]

    rad_angle = np.arctan(abs(x_diff/y_diff))

    # Must check the sign of x_diff and y_diff to add the proper angle at the end of the arctan calculation
    # This is needed since position angle is measured clockwise from north and are only calculating the angle 
    if (x_diff < 0) and (y_diff < 0):
        tail_angle = np.rad2deg(rad_angle) + 180

    if (x_diff < 0) and (y_diff > 0):
        tail_angle = 360 - np.rad2deg(rad_angle)

    if (x_diff > 0) and (y_diff < 0):
        tail_angle = np.rad2deg(rad_angle) + 90

    if (x_diff > 0) and (y_diff > 0):
        tail_angle = np.rad2deg(rad_angle)

    return tail_angle

def find_distance_conversions(name):
    '''
    Takes in the name of a galaxy and attempts to determine a conversion factor between the angular
    separation between two points in the sky and the physical distance between them in kpc.  If the
    galaxy name is not recognized under NASA Extragalactic Database, will arbitrarily give the proper
    conversion value for a galaxy in the Coma cluster.  Note that this is determined using the redshift
    of the galaxy in question. 
    
    Parameters
    ----------
    name: str
        the name of the galaxy one wishes to find the tail angle of   

    Returns
    -------
    distance_conversion_kpc_arcsec: float
        the conversion factor between angular separation and physical distances for an object, in units
        of kpc per arcsec
    '''

    # Attempts to pull the redshift data of the galaxy from NED.  If not possible, it'll use the redshift
    # of the Coma cluster (which all the RPS galaxies I'm studying are a part of)
    try:
        galaxyinfo = Ned.query_object(name)
        redshift = galaxyinfo[0]["Redshift"]
    
    except:
        redshift = 0.0234

    # Calculates the recessional velocity of the galaxy and, using Hubble's Law, determines a conversion factor
    recessional_velocity = (2.9979 * 10**5) * redshift
    proper_dist = recessional_velocity / 70
    distance_conversion_kpc_arcsec = (proper_dist / 206265) * 1000

    return distance_conversion_kpc_arcsec
    
        


def reconvert_distances_to_pixels(distance, conversion_factor, wcs, disk_com, ns = True, ne = True):
    '''
    Gives the pixel value of an image point that's a certain physical distance away from the galaxy 
    disk's centroid

    Parameters
    ----------
    dist: int
        the distance away from the galaxy centroid that one wishes to find the pixel value of, in kpc
    conversion_factor: float
        the conversion factor between kpc and arcsecs
    wcs: WCS object
        the WCS object of the galaxy image being used for the tail angle measurement
    disk_com: int array
        the pixel coordinates of the galaxy disk's flux-weighted centroid
    ns: bool (Default: True)
        sets whether the function finds the position north/south or east/west of the disk center of mass
    ne: bool (Default: True)
        sets whether the function finds the position goes north/east or south/west of the disk center of mass
       
    Returns
    -------
    x_pixel_position: int
        the x position of the point a certain distance out from the disk's centroid
    y_pixel_position: int
        the x position of the point a certain distance out from the disk's centroid
    '''

    # Converts the disk centroid coordinates to world coordinates
    com_wcs = wcs.pixel_to_world(disk_com[0] - 1, disk_com[1] - 1)
    
    # Checks if should calculate offset in ra or dec, then checks if should add or subtract that offset
    # Once has the new ra and decs, converts back to pixel values
    if ns == True:
        disk_y_com_wcs = com_wcs.dec.to(u.deg).value
        y_dist_offset = SkyCoord(0, distance / conversion_factor, unit = u.arcsec).dec.value

        if ne == True:
            new_dec = disk_y_com_wcs + y_dist_offset
            new_sky_coord = SkyCoord(com_wcs.ra, new_dec, unit = u.deg)

        else:
            new_dec = disk_y_com_wcs - y_dist_offset
            new_sky_coord = SkyCoord(com_wcs.ra, new_dec, unit = u.deg)

        x_pixel_position, y_pixel_position = wcs.world_to_pixel(new_sky_coord)

    else:
        disk_x_com_wcs = com_wcs.ra.to(u.deg).value
        x_dist_offset = SkyCoord(distance / conversion_factor, 0, unit = u.arcsec).ra.value
        if ne == True:
            new_ra = disk_x_com_wcs + x_dist_offset
            new_sky_coord = SkyCoord(new_ra, com_wcs.dec, unit = u.deg)

        else:
            new_ra = disk_x_com_wcs - x_dist_offset
            new_sky_coord = SkyCoord(new_ra, com_wcs.dec, unit = u.deg)

        x_pixel_position, y_pixel_position = wcs.world_to_pixel(new_sky_coord)

    return x_pixel_position, y_pixel_position

def set_slider_bounds(data, conversion_factor, wcs, disk_com, ns = True):
    '''
    Sets appropiate slider bounds based on the image provided.

    Parameters
    ----------
    data: 2d array
        the galaxy image data
    conversion_factor: float
        the conversion factor between kpc and arcsecs
    wcs: WCS object
        the WCS object of the galaxy image being used for the tail angle measurement
    disk_com: int array
        the pixel coordinates of the galaxy disk's flux-weighted centroid
    ns: bool (Default: True)
        sets whether the function finds the position north/south or east/west of the disk center of mass
 
    Returns
    -------
    upper_slider_bound: float
        the upper bound on the distance out from the galaxy disk allowed
    '''
    #Converts upper and lower pixel values in given direction to wcs, then multiplies by the conversion factor
    if ns == True:
        upper_coord = wcs.pixel_to_world(0,data.shape[0] - 1)
        lower_coord = wcs.pixel_to_world(0, disk_com[1])
        angular_separation = upper_coord.separation(lower_coord)
        upper_slider_bound = angular_separation.to(u.arcsec).value  * conversion_factor
    
    else:
        upper_coord = wcs.pixel_to_world(data.shape[1] - 1, 0)
        lower_coord = wcs.pixel_to_world(disk_com[0], 0)
        angular_separation = upper_coord.separation(lower_coord)
        upper_slider_bound = angular_separation.to(u.arcsec).value  * conversion_factor

    return upper_slider_bound

def determine_background(data):
    '''
    Calculates the average background value of the image.
    
    Parameters
    ----------
    data: 2d array
        the galaxy image one wishes to find the background at.  

    Returns
    -------
    background_value: float
        the average background value in the image
    '''

    eighth_bound_horizontal = int((data.shape[1] / 8))
    eighth_bound_vertical = int((data.shape[0] / 8))

    # Takes small slices from the bottom corners of the image and averages them
    bottom_left_corner = data[0:eighth_bound_vertical - 1, 0:eighth_bound_horizontal - 1]
    bottom_right_corner = data[0:eighth_bound_vertical - 1, data.shape[0] - eighth_bound_horizontal - 1]

    background_value = (np.mean(bottom_left_corner) + np.mean(bottom_right_corner)) / 2

    return background_value

def one_dim_center_of_mass(values, rows = True):
    '''
    Calculates the flux-weighted center of mass along each row/column of an image.
    
    Parameters
    ----------
    values: 2d array
        the galaxy image one wishes to use for the center of mass calculation.  
    rows: bool (default: True)
        determines whether to calculate the center of masses along each row or column.  

    Returns
    -------
    coms: float array
        the flux-weighted center of mass along each row/column
    '''
    coms = []
    vals = values

    if rows == False:
        vals = values.T

    for row in vals:
        norm = np.sum(row)
        weighted_vals = []

        for i in range(len(row)):
            weighted_pos = i * row[i]
            weighted_vals.append(weighted_pos)

        row_center_of_mass = np.sum(weighted_vals) / norm
        coms.append(row_center_of_mass)

    return coms

def two_dim_center_of_mass(values):
    '''
    Calculates the flux-weighted center of mass in a given 2-dimensional region.
    
    Parameters
    ----------
    values: 2d array
        the image data/region one wishes to find the center of mass of.  
    
    Returns
    -------
    com_coords: float array
        the exact flux-weighted center of mass coordinates
    '''
    data = [values.T, values]
    com_coords = []
    for vals in data:
        row_sums = []
        weighted_sums = []

        # Sums all the values along the row, then finds the center of mass of the
        # summed rows to get the y coordinate of the centroid (and repeats for the 
        # columns to get the x coordinate) 
        for row in vals:
            row_sum = np.sum(row)
            row_sums.append(row_sum)

        norm = np.sum(row_sums)

        for i in range(len(row_sums)):
            weighted_sum = i * row_sums[i]
            weighted_sums.append(weighted_sum)

        com = np.sum(weighted_sums) / norm
        com_coords.append(com)

    return com_coords   