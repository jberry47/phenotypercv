For the most up-to-date version of this program, please visit [DDPSC PhenotyperCV](https://github.com/jberry47/ddpsc_phenotypercv)

What is contained in this repo, is the exact code that was used to produce data found in the publication https://peerj.com/articles/5727/ 


# phenotypercv.cpp

Title: PhenotyperCV
Author: Jeffrey Berry (Bart Lab)
Description: This program is for segmenting and measuring plants from the Bellweather Phenotyping
Facility. Segmentation is achieved by supplying a background image that does not contain a plant
and using the difference between that and a supplied image to threshold on. Further processing is
done to remove artifacts that arise. After segmentation is complete, shapes and color profile are
reported in corresponding user-specified files.

Usage: There are three modes of use (VIS, VIS_CH, and NIR) and they correspond to what type of image
and processing that the program will perform. Depending on what is chosen, the required inputs change
1) VIS - Segment and measure plant in RGB images
   Example: ./PhenotyperCV VIS input_image.png background_image.png shapes.txt color.txt
2) VIS_CH - standardize, segment, and measure plant in RGB iamges (requires additional input, see below)
   Example: ./PhenotyperCV VIS_CH input_image.png background_image.png shapes.txt color.txt
NOTE Processing using the VIS_CH mode requries two additional items: a card_masks/ folder that contains
masks for each of the chips and target_homography.csv file that is the desired color space. The csv
file can be created for you using the provided tool, get_target_homography.cpp and redirecting the output,
but can also be manually obtained using your favorite image editor.
3) NIR - segment and measure plant in near-infrared images
   Example: ./PhenotyperCV NIR input_image.png nir_color.txt

Compiling Notes:
-I/usr/local/include/opencv -I/usr/local/include/opencv2 -I/usr/include/Eigen
-L/usr/local/lib lopencv_core -lopencv_features2d -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs

On Danforth Center Bioinformatics Infrastructure:
g++ -I/shares/bioinfo/installs/opencv-3.1.0/install/include -I/shares/bioinfo/installs/eigen/Eigen -L/shares/bioinfo/installs/opencv-3.1.0/install/lib -lopencv_imgproc -lopencv_imgcodecs -lopencv_core phenotypercv.cpp -o PhenotyperCV

To use the program in a parallelized fashion, pipe find into xargs with the -P flag followed by
number of cores like this:
   Usage VIS or VIS_CH:
      find Images/ -name 'VIS_SV*' | xargs -P8 -I{} ./PhenotyperCV VIS_CH {} background_image.png shapes.txt color.txt'
   
   Usage NIR:
      find Images/ -name 'NIR_SV*' | xargs -P8 -I{} ./PhenotyperCV NIR {} nir_color.txt



# avg_images.cpp

Title: Avg_Images
Author: Jeffrey Berry (Bart Lab)

Description:
Program takes single column stdout input of file paths to images to be averaged

Compiling notes:
same as phenotypercv.cpp except this does not require Eigen

Usage:
cat Images/SnapshotInfo.csv | grep Fm000Z | grep VIS_SV | awk -F'[;,]' '{print "Images/snapshot"$2"/"$12".png"}' | ./Avg_Images



# camera_calib.cpp

Title: Camera_Calib
Author: Jeffrey Berry (Bart Lab)

Description: 
This program extracts the RGB information of the color checker in the image. This requires a folder called card_masks that contains a binary image of each of the chips in the checker. See example folder. The output is written to stdout so you have to redirect it to target_homography.csv to be used by phenotypercv.cpp. 

Compiling notes:
same as phenotypercv.cpp except this does not require Eigen

Usage: 
./Camera_Calib average_images.png > target_homography.csv



# make_rois.txt

Title: Make_ROIs
Author: Jeffrey Berry (Bart Lab)

Description: This is an ImageJ macro that loops through ROI's in the ROI manager and creates the binary images required for camera_calib.cpp and phenotypercv.cpp (VIS_CH mode). 

Usage:
First you must manually make a ROI with the selector tool of your choosing and go through the chips one by one and add them to the ROI manager. The easiest way to do this is to open the manager first then move the ROI over a chip and press "t". Then move the next chip and again press "t" and so on. After all the ROIs are in the manager, go to Plugins > Macros > Run and execute this file. 
