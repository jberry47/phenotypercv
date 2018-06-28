#include <opencv2/opencv.hpp>
#include <cv.h>
#include <highgui.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

using namespace cv;
using namespace std;

float get_color(Mat img,Mat &mask){
	Mat img1 = img;
    Mat hist;
	int dims = 1;
	int histSize = 255;
	float hranges[] = { 0, 255 };
	const float *ranges = {hranges};

	calcHist(&img1,1,0,mask,hist, dims, &histSize, &ranges ,true ,false);

	int sum=0;
	for(int i = 0;i<256;i++){
		sum += hist.at<float>(i,0);
	}
	Mat weights = hist/sum;
	float hist_avg=0.0;
	for(int i = 0;i<256;i++){
		hist_avg += i*weights.at<float>(i,0);
	}
	return hist_avg;
}

int main(int argc, char** argv){
	Mat cc = imread(argv[1]);
    vector<Mat> bgr;
    split(cc, bgr);
    Mat b = bgr[0];
    Mat g = bgr[1];
    Mat r = bgr[2];
   	for(unsigned int i=1;i<23;i++){
   	    	stringstream ss;
   	    	ss << i;
   	    	string str = ss.str();
   	    	string file_name = "card_masks/"+str+"_mask.png";
   	    	Mat mask = imread(file_name,0);
   	    	Mat cc;
   	   		threshold(mask,cc,90,255,THRESH_BINARY);

   	    	float b_avg = get_color(b, cc);
   	    	float g_avg = get_color(g, cc);
   	    	float r_avg = get_color(r, cc);

   	    	//-- Write histogram averages to cout
   	    	cout  << r_avg << ","<< g_avg << "," << b_avg << endl;
   	}
}
