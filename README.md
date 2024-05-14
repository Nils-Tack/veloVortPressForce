# veloVortPressForce
A suite of master scripts used to plot velocity and vorticity fields using masks when an object is in the field of view. Experimentally derived pressure fields are also computed and used to measure the force components (thrust, drag, lateral component) along the object's outline.

The following sets of MATLAB scripts compute and plot velocity, vorticity and pressure fields and measure experimentally derived forces along an object in the fluid. Specifically, the scripts compute the instantaneous masks and outlines of the object/animal, as well as its centerlines (if it is an elongated object) to compute reliable pressure estimates around the body and measure the resulting forces. 

The scripts were written to be modular and can be reorganized as needed. Additionally, the scripts can be modified to perform calculations derived from the main features presented here. For instance, the centerlines are currently used to measure the heading of the organism in the field of view. However, the user can further develop the script by also finding the center of mass to perform torque calculations and determine the swimming speed of the object.

To perform the desired computations, proceed as follows:
1)	Create a folder containing the following subfolders:
-	Figures
-	Forces
-	Images
-	masksBW
-	outlines
-	PIVclean
-	PIVoriginal
-	pressure
-	centerlines

2)	Load the raw PIV data (as .txt files) in the ‘PIVoriginal folder and the corresponding images in the ‘images’ folder. Make sure the number of files is the same in these two folders. Note: When compouting velocity fields in DaVis, it is important to remember Davis assigns a velocity field on the first image in a sequence but cannot compute the velocity field from the last images. Warning: Always set the (0,0) origin of the axes in DaVis to the bottom left corner of the frame. To do so, DaVis has a setting in the scaling tab to align the y-axis origin. Because image axes origins are set to the top left corner by default, the user must offset the origin in DaVis with the height of the image in pixel.

3)	Load the ‘masterComputePressure’ script found in the ‘Master scripts’ folder.
-	Load the paths and make the directories.
-	Convert the raw PIV data to the .csv format for Queen2.
-	Generate and export the image masks to isolate the object. Select the desired options depending on your case. Also change some of the contrast and binarization thresholding parameters in the ‘makeMask3’ function if necessary
-	Load the scale (opts.scale) in px/mm (it will be converted to px/m) to scale the outlines properly.
-	Generate and export the outlines from the BW masks.
-	Repair the outlines manually if necessary. Follow the instructions at the top of the pop-up figure that appears upon initiating the ‘repairOutlines’ function.
-	Optional: export clean BW masks using the cleaned outlines.
-	Open the Queen2 ‘parameter2’ parameter file and change the parameters according to your case. Don’t forget to set your PIV and pressure file directories.
-	Run ‘Queen2’
-	Automatically move the pressure files Queen2 outputs from the ‘PIVclean’ folder to the ‘pressure’ folder.

4)	Load the ‘masterPlotVeloVortPress’ script found in the ‘Master scripts’ folder to plot the velocity, vorticity, and pressure fields.
-	Set the directories and paths.
-	Load the scale to scale the figure axes properly.
-	Import the velocity fields. This step takes a few seconds.
-	Convert sections of the velocity fields falling outside the computation domain to NaN. This would be the case if your version of DaVis outputs vectors = 0 outside of the geometric masks defining the region of interest for computing (if applicable).
-	Plot the mean velocity field and subtract the bulk flow if necessary. Note: change the opts.flowDirection if you need a different orientation or if there is no net bulk flow (1 = u (x direction); 2 = v (y direction); 3 = no flow).
-	Plot the velocity and vorticity fields. Run a test frame first to check that the options you selected match your expectations. Once satisfied, export all the frames. All the frames will be output in the ‘figures’ folder. The user can choose between exporting simple .jpg figures and/or .svg vector files that can be reworked in Illustrator. Default options: velocity vectors overlaid on top of the vorticity fields. The plots use the real image of the object, and the vorticity and velocity fields are masked using the outlines. These options can be changed to plot the mask instead of the image or plot only the vorticity fields without the vectors. Also change the plotting parameters to control the vorticity magnitude of the colorbar (symmetrical) and the density of vectors in the velocity fields.
-	Import the velocity fields. This step takes a few seconds.
-	Plot the pressure fields. The process is the same as for the velocity fields. Select the plotting parameters you need and run a test frame. Once satisfied, export all the pressure fields. All the frames will be output in the ‘figures’ folder. The user can choose between exporting simple .jpg figures and/or .svg vector files that can be reworked in Illustrator.

5)	Load the ‘masterCenterlineHeading’ script found in the ‘Master scripts’ folder to calculate the centerline and measure the heading of the object.
-	Set the directories and paths.
-	Load the scale to scale the figure axes properly.
-	Load the outlines or clean masksBW. If selecting the option to load the outlines, then a sub-function converts the outlines to BW images to perform morphological operations to skeletonize the masks.
-	Compute the centerlines. First, test a centerline to confirm that the centerlines are plotted properly. Then the script computes the centerlines for all the frames. At this stage, the centerline is re-scaled in pixels for convenience. Also, they may not be sorted properly, with the anterior point of the object identified as the first element in the centerline matrix.
-	Finalize the centerlines by identifying the anterior point as the first point of the centerline matrix. This stage of the computation also measures the heading of the object in the field of view. The heading is the swimming direction relative to the true north in the image (+y-axis, top edge of the frame). This function triggers a pop-up figure where the user will be asked to click on the anterior section of the object. A search area is thus defined around that point to automatically track the anterior point of the object in every frame. This local search finds the element corresponding to this anterior point to re-organize the centerline with that point as the first element of the centerline matrix. The search is assisted by predicting the displacement of the anterior point over time based on the displacement of the previous frames. In the event that the script cannot find the anterior point, the user will be prompted to select it manually to re-start the search. Options: the radius of the search area is set to 100 pixels. The size of the search area depends on the image resolution and the time increment between the frames, The quality of the outlines, and this thew resulting centerlines can affect the precision on the prediction to find the anterior point. This parameter can be changed in the ‘cleanCenterlines’ function found in the ‘Functions’ folder.
-	Convert the centerlines to meter.
-	Export the centerlines as individual .csv files for each frame in the ‘centerlines’ folder. This step also saves a matrix containing the heading of the object.
-	A last section allows you to plot individual centerlines along the mask that shows the identification of the anterior point and the position of the center of mass. In this case, the center of mass is located at a fixed distance from the anterior section, relative to the total length of the object. The default value is 30% of the length of the object.

6)	Load the ‘masterForcePlot’ script found in the ‘Master scripts’ folder to calculate and plot the force vectors. This function uses the heading of the object, provided either as a fixed, predetermined singular value or as a n×1 matrix containing the heading of the object relative to true north in degrees for each frame.
-	Set the directories and paths.
-	Loas all the options, including the scale, the time interval between frames, the heading option (fixed or dynamic). Swimming is currently considered fixed and is a legacy item that was used to compute power and the Froude efficiency. This can be changed in future iterations to accept dynamic changes in swimming speed based on the vector component of the displacement of the center of mass relative to the heading from one frame to the next. Other options included plotting preferences such as keeping the original image or simply plotting a clean mask to make a clean animation. The size of the vectors can also be changed by the factor opts.forceScale. If the reference scale is too long, it can be changed on line 12 of the ‘forceFromPressure’ function found in the ‘functions\forceComputation’ folder. This will also automatically update the reference vector label.
-	Test one frame and confirm that the options match your expectations. The function can accept four types of plots: 1) only plotting the axial, thrust and drag forces (relative to the heading), 2) only plotting the lateral forces, 3) Plotting the true force vector resulting from the lateral and axial components, and 4) A subplots showing cases 1–3. The figures can be saved as .jpg and/or .svg. The legend is not provided on to prevent cluttering the figure. However, it is plotted as a separate figure upon plotting a test frame.
