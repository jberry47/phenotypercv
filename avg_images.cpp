/* Title: Average list of images
 * Author: Jeffrey Berry (Bart Lab)
 *
 * Description:
 * Program takes single column stdout input of file paths to images to be averaged
 *
 * Usage:
 * cat SnapshotInfo.csv | grep Cc00AB | awk -F",|;" '{print "Images/snapshot"$2"/"$12".png"; print "Images/snapshot"$2"/"$13".png"}' | ~/cuda-workspace/Avg_Images/Debug/Avg_Images
 */

#include <opencv2/opencv.hpp>
#include <cv.h>

using namespace cv;
using namespace std;

int main(){
	string line;
	Mat avg;
	vector<Mat> avg_bgr(3);
	int counter = 0;
	while(cin) {
		if(getline(cin,line)) {
			if(counter == 0){
	    		avg=imread(line);
	    		avg.convertTo(avg, CV_64F);
    			split(avg,avg_bgr);
	    		counter++;
	    	}else{
	        	Mat inputImage = imread(line);
	        	inputImage.convertTo(inputImage, CV_64F);
	    		vector<Mat> in_bgr(3);

    			split(inputImage,in_bgr);
    			avg_bgr[0] = (avg_bgr[0]+in_bgr[0]);
    			avg_bgr[1] = (avg_bgr[1]+in_bgr[1]);
    			avg_bgr[2] = (avg_bgr[2]+in_bgr[2]);
	        	counter++;
	    	}
	    }
	}
	avg_bgr[0] = (avg_bgr[0])/counter;
	avg_bgr[1] = (avg_bgr[1])/counter;
	avg_bgr[2] = (avg_bgr[2])/counter;
	Mat adjImage;
	merge(avg_bgr,adjImage);
	adjImage.convertTo(adjImage, CV_64F);
	imwrite("average_images.png",adjImage);

	return 0;
}
