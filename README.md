# carps
Calculating the Angle of Ram Pressure Stripping! 

## Package Setup/Installation
If you wish to use CARP's web UI or any of CARP's functions for personal use, simply download the carps repository and run the following command in your command terminal.

```
pip install carps
```

Doing so will install the package onto your computer, given that you have the necessary dependencies for this package (which include numpy, matplotlib, astropy, astroquery, pandas, streamlit, and streamlit_image_coordinates).

From there, one can either locally run the web UI of this program using the instructions below or use any of the modules created for this program (such as data_analysis) by importing it with the following command:

```
import CARPS
```
## Website Setup 
Once the github repository is downloaded and pip installed, navigate to the main CARPS folder via the command termimal.

Once you are in the proper directory (which should contain the files "carps_site.py" and "data_analysis.py"), you can locally open the website by running the following command.

```
streamlit run carps_site.py
```

When this is done, your terminal should automatically open your preferred internet browser with the web ui all set up.  If not, check your command terminal for a url and copy and paste it into your web browser to open the website.


## How to use CARPS
1. Once the web interface is open, nothing particularly exciting will happen until you upload a FITS formatted file of the galaxy you'd like to measure the tail position angle of.  Note that this FITS file must contain WCS information in order for the program to work.

<p align="center">
<img src="https://github.com/rmunoz2002/carps/blob/main/readme_images/file_upload.PNG" width="500">
   </p>

2. Upon uploading a valid FITS file, you will be prompted to input the galaxy's name.  For the best experience, the inputted name should be any name for the galaxy recognized by the NASA Extragalactic Database.  If a non-recognized name is inputted, it will affect the program's ability to precisely convert WCS coordinates to physical distances (as it will default to a conversion factor for the center of the Coma Cluster).  

<p align="center">
<img src="https://github.com/rmunoz2002/carps/blob/main/readme_images/nameinput.PNG" width="500">
   </p>

3. Once a name is input it will display your FITS image (if you need to invert your image so high intensity corresponds to high pixel values, one can do so by clicking the "data correction" select box) and provide you with both a radio button widget and a drop down menu for pre-calculation image processing.  Begin by first selecting the radio button that correctly indicates the general direction the tail is pointing in (north, south, east, or west).  If the galaxy's tail is pointing in between these cardinal directions, one can simply select either of the cardinal directions associated with the tail based on whether its pointing more in one direction or the other.

<p align="center">
<img src="https://github.com/rmunoz2002/carps/blob/main/readme_images/radiowidget.PNG" width="500">
   </p>

4. From there, one must define the region of the image associated with the actual galaxy's disk.  This is done by simply clicking on the displayed image to place a point until you have completely bounded the galaxy disk.  If you misplace a point, you can simply remove it by pressing the "Remove Last Position" button.  Once the disk is bound to your liking, simply press the "Confirm Head" button.

<p align="center">
<img src="https://github.com/rmunoz2002/carps/blob/main/readme_images/diskdefine.PNG" width="500">
   </p>

5. This should add a red star to the image, which indicates the flux-weighted centroid of the galaxy itself.  From there, one can optionally continue image processing by changing the drop down menu to "Crop" and selecting areas that you wish to remove from the image (which will help with the final calculation's accuracy).  Note that one can select points outside of the bounds of the galaxy image itself and the program will correct those points when cropping.

<p align="center">
<img src="https://github.com/rmunoz2002/carps/blob/main/readme_images/cropdef.PNG" width="500">
   </p>

6. If the disk is defined, one can then select the distance out from the disk centroid to measure the tail angle at.  Note that these bounds on the distance are determined by the redshift of the galaxy and are thus subject to disagreement.  Regardless, once distance selection is done, pressing the "Calculate Tail Angle!" button will both update the graph to show the tail line (along with tail lines corresponding to just above and below the measured distance) and create a table with information surrounding the tail's position angle and its associated error.

<p align="center">
<img src="https://github.com/rmunoz2002/carps/blob/main/readme_images/selectdist.PNG" width="500">
   </p>

7. Once the table is filled out to your liking, press the "Export Data!" button to export the table as a .csv file. 

<p align="center">
<img src="https://github.com/rmunoz2002/carps/blob/main/readme_images/carpstable.PNG" width="500">
   </p>

<p align="center">
<img src="https://github.com/rmunoz2002/carps/blob/main/readme_images/finalimagecarps.PNG" width="500">
   </p>

Have fun!

