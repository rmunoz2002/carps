import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
from astropy.io import fits
from streamlit_image_coordinates import streamlit_image_coordinates
from PIL import Image, ImageDraw
from astropy.wcs import WCS
import math
import pandas as pd
from data_analysis import *

st.title('CARPS: Calculate the Angle of Ram Pressure Stripping!')

uploaded_file = st.file_uploader('RPS Galaxy File', type=['fit','fits','FIT','FITS'], accept_multiple_files=False, key='fileupload', help='Science data of a RPS galaxy you wish to find the tail angle of.  Accepts standard FITS formatting.', on_change=None, args=None, kwargs=None, disabled=False, label_visibility="visible")

janky_data = st.checkbox("Data Correction", value=False, key=None, help="Check this if your data is inverted or otherwise formatted in a janky/unorthodox way.", on_change=None, args=None, kwargs=None, disabled=False, label_visibility="visible")

# Initialize session_state variables if not already
if 'upperbound' not in st.session_state:
    st.session_state.upperbound = None

if 'lowerbound' not in st.session_state:
    st.session_state.lowerbound = None

if 'tail_coord' not in st.session_state:
    st.session_state['tail_coord'] = None

if 'data' not in st.session_state:
    st.session_state['data'] = None

if 'disk_com' not in st.session_state:
    st.session_state['disk_com'] = None

if 'background' not in st.session_state:
    st.session_state["background"] = None

if 'distance' not in st.session_state:
    st.session_state["distance"] = None

if 'upper_tail_angle' not in st.session_state:
    st.session_state["upper_tail_angle"] = []

if 'lower_tail_angle' not in st.session_state:
    st.session_state["lower_tail_angle"] = []

if 'tail_angle' not in st.session_state:
    st.session_state["tail_angle"] = []

if "distance_list" not in st.session_state:
    st.session_state["distance_list"] = []

if 'points' not in st.session_state:
    st.session_state["points"] = []

# FILE UPLOADING CHECKS
if uploaded_file is not None:
    # Checks to see if edits have already been made to the data (i.e. crops), and if so, uses edited data instead.
    if st.session_state["data"] is not None:
        galaxy_data = st.session_state["data"]
        hdu = fits.open(uploaded_file)
        wcs = WCS(hdu[0].header) 
    else:
        #Attempts to use the uploaded data file if no pre-existing data, but will bail if doesn't have WCS info 
        try:
            hdu = fits.open(uploaded_file)
            galaxy_data = hdu[0].data
            wcs = WCS(hdu[0].header)

            if janky_data:
                galaxy_data = -1 * galaxy_data + np.max(galaxy_data) + 1
        
        except:
            st.warning("File must have WCS info!")
            st.stop()

# If no file has been uploaded, does not continue with the code
else:
    st.warning("Please input a file! xoxo")
    st.stop()


nameinput = st.text_input("Galaxy Name", value="", max_chars=None, key=None, type="default", help="Name of the galaxy. Name needs to be a valid name in the NASA Extragalactic Database", autocomplete=None, on_change=None, args=None, kwargs=None, placeholder=None, disabled=False, label_visibility="visible")

# NAME UPLOADING CHECK
if nameinput:
     galaxy_name = nameinput

else:
    st.warning("Please input a name tee hee")
    st.stop()


# PLOTTING FUNCTIONS
def remove_last_position():
    '''
    Simply removes the last added point from the state session.
    '''
    if len(st.session_state["points"]) == 0:
        st.write("No points to remove!") 

    else:
        st.session_state["points"].pop()

def get_ellipse_coords(point: tuple[int, int]) -> tuple[int, int, int, int]:
    '''
    Gets the coordinates of a circular point centered about a given coordinate. 
    
    Parameters
    ----------
    point: int tuple
        the coordinates of the circle center.  
    
    Returns
    -------
    ellipse_coords: 4-tuple
        the pixel coordinates of vertices and covertices of the circle.
    '''
    center = point
    radius = 3
    return (
        center[0] - radius,
        center[1] - radius,
        center[0] + radius,
        center[1] + radius,
    )

def define_galaxy_disk(positions, imagefile, data):
    '''
    Defines a given region bounded within the coordinates selected with streamlit_image_coordinates. 
    
    Parameters
    ----------
    positions: int array
        the coordinates of the vertices that will bound the region.  
    imagefile: str
        the path to the image file being used to display the image on the website.
    data: 2d array
        the actual image data being used for the tail angle calculation
    
    Returns
    -------
    galaxy_disk_region: 2d array (bool)
        the actual data mask of the region bounded by the inputted points
    '''

    # Checks if there's sufficient points to bound a 2D region
    if len(positions) < 3:
        st.write("Not enough points selected!")
        
    else:
        # Takes the image file, converts it to a value array, then finds its shape along with the 
        # shape of the actual data 
        with Image.open(imagefile) as input:
            newpoints = []
            inputdata = np.asarray(input)

            imageshape = inputdata.shape
            datashape = data.shape

            # Calculates the difference in shapes between the two and finds the necessary correction 
            x_positional_offset = int((imageshape[1] - datashape[1]) / 2)
            y_positional_offset = int((imageshape[0] - datashape[0]) / 2)

            # Corrects each inputted position to correspond to its proper pixel value in the data array
            for position in positions:
                raw_x_position = position[0] - x_positional_offset
                raw_y_position = imageshape[0] - position[1] - y_positional_offset
                new_point = raw_x_position, raw_y_position
                newpoints.append(new_point)

            # Checks if a point is outside the bounds of the actual data array and, if so, snaps its
            # out of bound position(s) to the nearest position that is in bounds
            for new_coord in newpoints:
                if  new_coord[0] < 0:
                    y = list(new_coord)
                    y[0] = 0
                    new_coord = tuple(y)

                if new_coord[0] > datashape[1] - 1:
                    y = list(new_coord)
                    y[0] = datashape[1] - 1
                    new_coord = tuple(y)

                if  new_coord[1] < 0:
                    y = list(new_coord)
                    y[1] = 0
                    new_coord = tuple(y)

                if new_coord[1] > datashape[0] - 1:
                    y = list(new_coord)
                    y[1] = datashape[0] - 1
                    new_coord = tuple(y)
            
            # Creates a polygonal path using these 'corrected' points along with a meshgrid of all the coordinate
            # combinations/positions. 
            mask_path = mpath.Path(newpoints)

            x = np.linspace(0,datashape[1] - 1, datashape[1])
            y = np.linspace(0,datashape[0] - 1, datashape[0])
            xv, yv = np.meshgrid(x, y)

            coordinates = np.vstack([xv.ravel(), yv.ravel()]).T

            # Finally creates a mask of all the values within the polygonal path.
            galaxy_disk_region = np.reshape(mask_path.contains_points(coordinates),(datashape[0],datashape[1]))

            return galaxy_disk_region
               
# If a file has been successfully uploaded and converted into a data array, display the image and 
# allow for user to continue onto next steps in the tail-angle measurement process.
if galaxy_data is not None:
    editmode = st.selectbox(label='Edit Mode', options = ['Select Galaxy Disk (Required)','Crop (Optional)'], index=0, key=None, help=None, on_change=None, args=None, kwargs=None, disabled=False, label_visibility="visible")
    fig = plt.figure()
    fig.add_subplot(projection=wcs) 

    # Make the graph appropiately color mapped. 
    if st.session_state['background'] is not None:
        plt.imshow(galaxy_data, cmap = 'gray', vmin = st.session_state['background'], vmax = np.max(galaxy_data))

    else:
        plt.imshow(galaxy_data, cmap = 'gray', vmin = np.min(galaxy_data), vmax = np.max(galaxy_data))

    # If a disk centroid has been found, display it.
    if st.session_state['disk_com'] is not None:
        plt.scatter(st.session_state['disk_com'][0], st.session_state['disk_com'][1], marker = '*', color = 'r', label = 'Galaxy Disk Centroid')
        plt.legend()

    # If tail coordinates (and in turn tail angles) have been calculated, display those points 
    # and determine error bounds on the tail angle.
    if st.session_state['tail_coord'] is not None:
        plt.plot([st.session_state['disk_com'][0], st.session_state['tail_coord'][0]], [st.session_state['disk_com'][1], st.session_state['tail_coord'][1]], label = 'Galaxy Tail Line')
        plt.plot([st.session_state['disk_com'][0], st.session_state['upperbound'][0]], [st.session_state['disk_com'][1], st.session_state['upperbound'][1]], label = 'Outer Bound')
        plt.plot([st.session_state['disk_com'][0], st.session_state['lowerbound'][0]],[st.session_state['disk_com'][1], st.session_state['lowerbound'][1]], label = 'Inner Bound')
        plt.legend(fontsize="7")

        
        stds = []

        # Standard deviation calculation for each tail angle measured (added in quadrature).
        for i in range(len(st.session_state.tail_angle)):
            stds.append(((st.session_state["lower_tail_angle"][i] - st.session_state["tail_angle"][i])**2 + (st.session_state["upper_tail_angle"][i] - st.session_state["tail_angle"][i])**2)**(1/2))

        # Make a Pandas dataframe containing all the appropiate information for each measurement and display it.
        tail_info = {'Galaxy Name': galaxy_name,'Distance from Galaxy Disk (kpc)': st.session_state.distance_list, 'Position Angle of Tail (deg)': st.session_state["tail_angle"], 'Inner Bound on Tail Angle (deg)': st.session_state["lower_tail_angle"], 'Outer Bound on Tail Angle (deg)': st.session_state["upper_tail_angle"], 'Error on Tail Angle (deg)': stds}
        tail_dataframe = pd.DataFrame(data = tail_info)
        
        st.write(tail_dataframe)

        #Allow for users to export their data as a csv file
        if st.button('Export Data!', key=None, help=None, on_click=None, args=None, kwargs=None, type="secondary", disabled=False, use_container_width=False):
            tail_dataframe.to_csv('tail_data.csv')

    tail_dir = st.radio("General Tail Direction", ("North", "South", "East", "West"), index=0, key=None, help=None, on_change=None, args=None, kwargs=None, disabled=False, horizontal=False, label_visibility="visible")
    
    # Check for tail direction for future calculations 
    if tail_dir == "North":
        ns = True
        ne = True
        
    if tail_dir == "South":
        ns = True
        ne = False

    if tail_dir == "East":
        ns = False
        ne = True

    if tail_dir == "West":
        ns = False
        ne = False
    
    # Update the image being displayed.
    plt.savefig('imagedata.jpg', overwrite = 'True') 

    if editmode == 'Select Galaxy Disk (Required)':
        with Image.open("imagedata.jpg") as img:
            draw = ImageDraw.Draw(img)
            
            if st.button(label = 'Remove Last Point'):
                remove_last_position()
                st.experimental_rerun()

            if st.button(label = 'Confirm Head', key=None, help=None, on_click=None, args=None, kwargs=None, type="secondary", disabled=False, use_container_width=False):
                
                disk_mask = define_galaxy_disk(st.session_state["points"],"imagedata.jpg", galaxy_data)
                galaxy_head = galaxy_data

                # Mask out everything but the galaxy disk itself
                inverse_disk_mask = np.logical_not(disk_mask)
                galaxy_head[inverse_disk_mask] = 1

                # Find the disk centroid
                disk_center_of_mass = two_dim_center_of_mass(galaxy_head)
                
                # Reset the selected points and update the disk centroid coordinates
                st.session_state['disk_com'] = disk_center_of_mass
                st.session_state['points'] = []
                st.experimental_rerun()
                    


            if len(st.session_state['points']) >= 3:
                draw.polygon(st.session_state['points'], outline= 'red')
            
            for point in st.session_state["points"]:
                coords = get_ellipse_coords(point)
                draw.ellipse(coords, fill="red")

            value = streamlit_image_coordinates(img, key="pil")

            if value is not None:
             point = value["x"], value["y"]

             if point not in st.session_state["points"]:
                    st.session_state["points"].append(point)
                    st.experimental_rerun()

    if editmode == 'Crop (Optional)':
        with Image.open("imagedata.jpg") as img:
            draw = ImageDraw.Draw(img)

            if st.button(label = 'Remove Last Point'):
                remove_last_position()
                st.experimental_rerun()

            if st.button(label = 'Crop Selected Area', key=None, help=None, on_click=None, args=None, kwargs=None, type="secondary", disabled=False, use_container_width=False):
                # Determine the region to crop out and find the average background value.
                crop_mask = define_galaxy_disk(st.session_state["points"],"imagedata.jpg", galaxy_data)
                image_background = determine_background(galaxy_data)
                
                # Set the cropped region to the background value, then update the data,
                # reset the points, and save the background value.
                galaxy_data[crop_mask] = image_background

                st.session_state["data"] = galaxy_data
                st.session_state["points"] = []
                st.session_state["background"] = image_background

                st.experimental_rerun()
                    
            # Display every point and, if applicable, draw a polygon defined by said points.
            if len(st.session_state['points']) >= 3:
                draw.polygon(st.session_state['points'], outline= 'red')
            
            for point in st.session_state["points"]:
                coords = get_ellipse_coords(point)
                draw.ellipse(coords, fill="red")

            # Getting coordinate values from user's clicks on the displayed image.
            value = streamlit_image_coordinates(img, key="pil")

            if value is not None:
             point = value["x"], value["y"]

             if point not in st.session_state["points"]:
                    st.session_state["points"].append(point)
                    st.experimental_rerun()


if galaxy_data is not None:
    if galaxy_name is not None:
        
        # If a disk centroid is defined, allow users to continue with the tail measurement process
        # by creating a distance slider with appropiate distance values.
        if st.session_state.disk_com is not None:
            disk_com = st.session_state.disk_com

            dist_conversion = find_distance_conversions(galaxy_name)
            upper_slider = set_slider_bounds(galaxy_data, dist_conversion, wcs, disk_com, ns)

            st.session_state.distance = st.slider('Select Distance from Galaxy Disk (kpc)', min_value = 1, max_value = math.floor(upper_slider) - 2, value = 0, step=1)

# Check to make sure everything necessary is there to calculate tail angle, then begin the tail angle
# calculation process.
if st.button(label = 'Calculate Tail Angle!'):
    if st.session_state["disk_com"] is None:
        st.write("You forgot to define the galaxy disk you silly goose :)")

    else:
        if st.session_state['distance'] is None:
            st.write("Gotta select a distance for that tail angle my dude.")

        else:
            # Get center of masses at each x or y coordinate and get the proper pixel coords a certain distance
            # out from the centroid and apply the proper index to the coms array to get the tail coordinate.
            coms = one_dim_center_of_mass(galaxy_data, ns)
            x_tail_coord, y_tail_coord = reconvert_distances_to_pixels(st.session_state.distance, dist_conversion, wcs, disk_com, ns, ne)
            
            # Get the proper tail coordinate based on general direction of the tail.
            if ns == True:
                y_index = int(y_tail_coord)
                tail_coordinate = [int(coms[y_index]), y_index]

            else:
                x_index = int(x_tail_coord)
                tail_coordinate = [x_index, int(coms[x_index])]

            st.session_state.tail_coord = tail_coordinate

            # Repeat the above but for slightly above and below the inputted distance for error estimation.
            upper_bound_x, upper_bound_y = reconvert_distances_to_pixels(st.session_state.distance + 2, dist_conversion, wcs, disk_com, ns, ne)
            lower_bound_x, lower_bound_y = reconvert_distances_to_pixels(st.session_state.distance - 2, dist_conversion, wcs, disk_com, ns, ne)

            if ns == True:
                upper_coord = [int(coms[int(upper_bound_y)]), int(upper_bound_y)]
                lower_coord = [int(coms[int(lower_bound_y)]), int(lower_bound_y)]
           
            else:
                upper_coord = [int(upper_bound_x), int(coms[int(upper_bound_x)])]
                lower_coord = [int(lower_bound_x), int(coms[int(lower_bound_x)])]
            
            distance = st.session_state.distance
            st.session_state["distance_list"].append(distance)

            # With the coordinates of the tail positions, find the corresponding position angles.
            tail_angle = find_angle_wrt_north(tail_coordinate, disk_com)
            tail_angle_upper = find_angle_wrt_north(upper_coord, disk_com)
            tail_angle_lower = find_angle_wrt_north(lower_coord, disk_com)

            # Update the session state by adding the calculated tail values, coords, and distance.
            st.session_state.upper_tail_angle.append(tail_angle_upper) 
            st.session_state.lower_tail_angle.append(tail_angle_lower)
            st.session_state.tail_angle.append(tail_angle)
            st.session_state.upperbound = upper_coord
            st.session_state.lowerbound = lower_coord
            st.experimental_rerun()
    
