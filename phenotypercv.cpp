/* Title: PhenotyperCV
 * Author: Jeffrey Berry (Bart Lab)
 *
 * Description: This program is for segmenting and measuring plants from the Bellweather Phenotyping
 * Facility. Segmentation is achieved by supplying a background image that does not contain a plant
 * and using the difference between that and a supplied image to threshold on. Further processing is
 * done to remove artifacts that arise. After segmentation is complete, shapes and color profile are
 * reported in corresponding user-specified files.
 *
 * Usage: There are three modes of use (VIS, VIS_CH, and NIR) and they correspond to what type of image
 * and processing that the program will perform. Depending on what is chosen, the required inputs change
 * 	1) VIS - Segment and measure plant in RGB images
 * 		Example: ./PhenotyperCV VIS input_image.png background_image.png shapes.txt color.txt
 *
 * 	2) VIS_CH - standardize, segment, and measure plant in RGB iamges (requires additional input, see below)
 * 		Example: ./PhenotyperCV VIS_CH input_image.png background_image.png shapes.txt color.txt
 * NOTE Processing using the VIS_CH mode requries two additional items: a card_masks/ folder that contains
 * masks for each of the chips and target_homography.csv file that is the desired color space. The csv
 * file can be created for you using the provided tool, get_target_homography.cpp and redirecting the output,
 * but can also be manually obtained using your favorite image editor.
 *
 * 	3) NIR - segment and measure plant in near-infrared images
 * 		Example: ./PhenotyperCV NIR input_image.png nir_color.txt
 * * g++ -I/shares/bioinfo/installs/opencv-3.1.0/install/include -I/shares/bioinfo/installs/eigen/Eigen -L/shares/bioinfo/installs/opencv-3.1.0/install/lib -lopencv_imgproc -lopencv_imgcodecs -lopencv_core phenotypercv.cpp -o PhenotyperCV
 *
 * Compiling Notes:
 * -I/usr/local/include/opencv -I/usr/local/include/opencv2 -I/usr/include/Eigen
 * -L/usr/local/lib lopencv_core -lopencv_features2d -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs
 *
 * On Danforth Center Bioinformatics Infrastructure:
 * g++ -I/shares/bioinfo/installs/opencv-3.1.0/install/include -I/shares/bioinfo/installs/eigen/Eigen -L/shares/bioinfo/installs/opencv-3.1.0/install/lib -lopencv_imgproc -lopencv_imgcodecs -lopencv_core phenotypercv.cpp -o PhenotyperCV
 *
 * To use the program in a parallelized fashion, pipe find into xargs with the -P flag followed by
 * number of cores like this:
 * 	Usage VIS or VIS_CH:
 * 		find Images/ -name 'VIS_SV*' | xargs -P8 -I{} ./PhenotyperCV VIS_CH {} background_image.png shapes.txt color.txt'
 * 	Usage NIR:
 * 		find Images/ -name 'NIR_SV*' | xargs -P8 -I{} ./PhenotyperCV NIR {} nir_color.txt
 *
 */

#include <opencv2/opencv.hpp>
#include <opencv/cv.h>
#include <highgui.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <opencv2/features2d.hpp>
#include <Eigen/Dense>


using namespace cv;
using namespace std;
using namespace Eigen;

int is_oof(Mat img){
	//-- Get contours of mask
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
    findContours( img, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE, Point(0, 0) );

    //-- Get contours of rectangular roi
    Mat src = Mat::zeros(img.size(),img.type())+255;

    vector<vector<Point> > contours_roi;
    vector<Vec4i> hierarchy_roi;
    findContours( src, contours_roi, hierarchy_roi, RETR_TREE, CHAIN_APPROX_SIMPLE, Point(0, 0) );

    int check = 0;
    //-- Keep only those contours that have a point inside roi
    for(unsigned int i=0; i < contours.size(); i++){
      	for(unsigned int j=0; j<contours[i].size(); j++){
      		int test = pointPolygonTest(contours_roi[0],Point2f(contours[i][j]),false);
      		if(test == 0){
      			check = 1;
      		}
       	}
    }
	return check;
}

vector<Point> keep_roi(Mat img,Point tl, Point br, Mat &mask){
	//-- Get contours of mask
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
    findContours( img, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE, Point(0, 0) );

    //-- Get contours of rectangular roi
    Mat src = Mat::zeros(img.size(),img.type());
    rectangle(src,tl,br,255,CV_FILLED);

    vector<vector<Point> > contours_roi;
    vector<Vec4i> hierarchy_roi;
    findContours( src, contours_roi, hierarchy_roi, RETR_TREE, CHAIN_APPROX_SIMPLE, Point(0, 0) );

    //-- Keep only those contours that have a point inside roi
    vector<Point> cc;
    Mat kept = Mat::zeros(img.size(),img.type());
    for(unsigned int i=0; i < contours.size(); i++){
      	for(unsigned int j=0; j<contours[i].size(); j++){
      		int test = pointPolygonTest(contours_roi[0],Point2f(contours[i][j]),false);
      		if(test==1 || test == 0){
      			for(unsigned int k=0; k<contours[i].size(); k++){
      				cc.push_back(contours[i][k]);
      			}
      			drawContours(kept, contours, i, 255, CV_FILLED);
      			break;
      		}
       	}
    }
	Mat kept_mask;
	bitwise_and(img,kept,kept_mask);

    mask = kept_mask;
	return cc;
}

float get_fd(Mat mask){
	//-- Need to remap the image to 512x512 so box counting can be used
	Mat img_bc;
	resize(mask, img_bc, Size(512,512), 0, 0, INTER_LINEAR);

	//-- Initializing variables
	double width = 512.0;
	double p = log(width)/log(double(2.0));
	VectorXf N = VectorXf::Zero(int(p)+1);
	double sumImg = sum(img_bc)[0];
	N(int(p)) = sumImg;

	//-- Boxcounting
	double siz;
	double siz2;
	float running_sum;
    for (int g = int(p)-1; g > 0; g--){
    	siz = pow(2.0, double(p-g));
    	siz2 = round(siz/2.0);
    	running_sum = 0;
    	for (int i = 0; i < int(width-siz+1); i = i+int(siz)){
    		for (int j = 0; j < int(width-siz+1); j = j+int(siz)){
    			img_bc.at<uchar>(i,j) = (bool(img_bc.at<uchar>(i,j)) || bool(img_bc.at<uchar>(i+siz2,j))
    				|| bool(img_bc.at<uchar>(i,j+siz2)) || bool(img_bc.at<uchar>(i+siz2,j+siz2)));
    			running_sum = running_sum+float(img_bc.at<uchar>(i,j));
    		}
    	}
    	N(g) = running_sum;
	}
    N = N.colwise().reverse().eval();

    //-- Getting bin sizes
    VectorXf R = VectorXf::Zero(int(p)+1);
    R(0) = 1.0;
    for (int k = 1; k < R.size(); k++){
    	R(k) = pow(2.0, double(k));
    }

    //-- Calculating log-log slopes
	float slope [R.size()-1];
	for(int i=1;i < R.size()-1 ;i++){
		slope[i] = (log10(N(i+1))-log10(N(i)))/(log10(R(i+1))-log10(R(i)));
	}

	//-- Getting average slope (fractal dimension)
	float sum = 0.0, average;
	int s_count =0;
	for(int i=1; i < R.size()-1; i++){
		if(-slope[i] < 2 && -slope[i] > 0){
			sum += -slope[i];
			s_count++;
		}
	}
	average = sum / s_count;
	return average;
}

vector<double> get_shapes(vector<Point> cc,Mat mask){
    //-- Get measurements
    Moments mom = moments(mask,true);
    double area = mom.m00;
    vector<Point>hull;
    convexHull( Mat(cc), hull, false );
    double hull_verticies = hull.size();
    double hull_area = contourArea(Mat(hull));
    double solidity = area/hull_area;
    double perimeter = arcLength(Mat(cc),false);
    double cmx = mom.m10 / mom.m00;
    double cmy = mom.m01 / mom.m00;
    Rect boundRect = boundingRect( cc );
    double width = boundRect.width;
    double height = boundRect.height;
    double circ = 4*M_PI*area/(perimeter*perimeter);
    double angle = -1;
    double ex = -1;
    double ey = -1;
    double emajor = -1;
    double eminor = -1;
    double eccen = -1;
    double round = -1;
    double ar = -1;
    if(cc.size() >= 6){
        Mat pointsf;
    	Mat(cc).convertTo(pointsf, CV_32F);
   	    RotatedRect ellipse = fitEllipse(pointsf);
   	    angle = ellipse.angle;
  	    ex = ellipse.center.x;
   	    ey = ellipse.center.y;
   	    if(ellipse.size.height > ellipse.size.width){
   	    	emajor = ellipse.size.height;
   	    	eminor = ellipse.size.width;
   	    }else{
   	    	eminor = ellipse.size.height;
  	   	    emajor = ellipse.size.width;
   	    }
   	    eccen = sqrt((1- eminor / emajor)*2);
   	    round = eminor/emajor;
   	    ar = emajor/eminor;
    }
    float fd = get_fd(mask);
    double oof = is_oof(mask);
    double shapes[20] = {area,hull_area,solidity,perimeter,width,height,cmx,cmy,hull_verticies,ex,ey,emajor,eminor,angle,eccen,circ,round,ar,fd,oof};
    vector<double> shapes_v(shapes,shapes+20);
    return shapes_v;
}

Mat get_color(Mat img,Mat mask){
	Mat composite;
	cvtColor(img,composite,COLOR_BGR2HSV);
    vector<Mat> channels1;
    split(composite, channels1);
    Mat hist;
	int dims = 1; // Only 1 channel, the hue channel
	int histSize = 180; // 180 bins, actual range is 0-360.
	float hranges[] = { 0, 180 }; // hue varies from 0 to 179, see cvtColor
	const float *ranges = {hranges};

	//-- Compute the histogram
	calcHist(&channels1[0],1,0,mask,hist, dims, &histSize, &ranges	,true ,false);
	return hist;
}

Mat get_nir(Mat img,Mat mask){
    Mat hist;
	int dims = 1;
	int histSize = 255;
	float hranges[] = { 0, 255 };
	const float *ranges = {hranges};

	//-- Compute the histogram
	calcHist(&img,1,0,mask,hist, dims, &histSize, &ranges	,true ,false);
	return hist;
}

float extractRGB_chips(Mat img,Mat &mask){
	//-- Averages the histogram for a given channel
	Mat img1 = img;
    Mat hist;
	int dims = 1;
	int histSize = 255;
	float hranges[] = { 0, 255 };
	const float *ranges = {hranges};

	calcHist(&img1,1,0,mask,hist, dims, &histSize, &ranges ,true ,false);

	int sum=0;
	for(int i = 0;i<255;i++){
		sum += hist.at<float>(i,0);
	}
	Mat weights = hist/sum;
	float hist_avg=0.0;
	for(int i = 0;i<255;i++){
		hist_avg += i*weights.at<float>(i,0);
	}
	return hist_avg;
}

MatrixXd getRGBarray(Mat img){
	//-- Loops over chips and gets RGB values of each one
	MatrixXd sourceColors(22,3);
	vector<Mat> bgr;
	split(img, bgr);
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
	   		float b_avg = extractRGB_chips(b, cc);
	   	    float g_avg = extractRGB_chips(g, cc);
	   	    float r_avg = extractRGB_chips(r, cc);
	   	    sourceColors(i-1,0) = b_avg;
	   	    sourceColors(i-1,1) = g_avg;
	   	    sourceColors(i-1,2) = r_avg;
	}
    return(sourceColors);
}

void get_standardizations(Mat img, float &det, MatrixXd &rh,MatrixXd &gh,MatrixXd &bh){
	//-- Extending source RGB chips to squared and cubic terms
	MatrixXd source1, source2, source3;
	source1 = getRGBarray(img);
	source2 = (source1.array() * source1.array()).matrix();
	source3 = (source2.array() * source1.array()).matrix();
	MatrixXd source(source1.rows(),source1.cols()+source2.cols()+source3.cols());
	source << source1, source2, source3;

	//-- Computing Moore-Penrose Inverse
	MatrixXd M = (source.transpose()*source).inverse()*source.transpose();

	//-- Reading target homography
	MatrixXd target(22,3);
	fstream file;
	file.open("target_homography.csv");
	string value;
	int rowCounter = 0;
	while ( getline ( file, value) )
	{
	     vector<float> result;
	     stringstream substr(value);
	     string item;
	     while (getline(substr, item, ',')) {
		     const char *cstr = item.c_str();
		     char* pend;
		     float num = strtof(cstr,&pend);
	         result.push_back(num);
	     }
	     target(rowCounter,0) = result[0];
	     target(rowCounter,1) = result[1];
	     target(rowCounter,2) = result[2];
	     rowCounter++;
	}

	//-- Computing linear target RGB standardizations
	rh = M*target.col(2);
	gh = M*target.col(1);
	bh = M*target.col(0);

	//-- Extending target RGB chips to squared and cubic terms
	MatrixXd target1, target2, target3;
	target2 = (target.array() * target.array()).matrix();
	target3 = (target2.array() * target.array()).matrix();

	//-- Computing square and cubic target RGB standardizations
	MatrixXd r2h,g2h,b2h,r3h,g3h,b3h;
	r2h = M*target2.col(2);
	g2h = M*target2.col(1);
	b2h = M*target2.col(0);
	r3h = M*target3.col(2);
	g3h = M*target3.col(1);
	b3h = M*target3.col(0);

	//-- Computing D
	MatrixXd H(9,9);
	H << rh.col(0),gh.col(0),bh.col(0),r2h.col(0),g2h.col(0),b2h.col(0),r3h.col(0),g3h.col(0),b3h.col(0);
	det = H.transpose().determinant();
}

Mat color_homography(Mat img, MatrixXd r_coef,MatrixXd g_coef,MatrixXd b_coef){
	Mat b, g, r, b2, g2, r2, b3, g3, r3;
	vector<Mat> bgr(3);
	split(img,bgr);

	//-- Computing linear, squared, and cubed images
	b = bgr[0];
	g = bgr[1];
	r = bgr[2];
	b2 = b.mul(b);
	g2 = g.mul(g);
	r2 = r.mul(r);
	b3 = b2.mul(b);
	g3 = g2.mul(g);
	r3 = r2.mul(r);

	//-- Computing homography
	b = 0+r*b_coef(0,0)+g*b_coef(1,0)+b*b_coef(2,0)+r2*b_coef(3,0)+g2*b_coef(4,0)+b2*b_coef(5,0)+r3*b_coef(6,0)+g3*b_coef(7,0)+b3*b_coef(8,0);
	g = 0+r*g_coef(0,0)+g*g_coef(1,0)+b*g_coef(2,0)+r2*g_coef(3,0)+g2*g_coef(4,0)+b2*g_coef(5,0)+r3*g_coef(6,0)+g3*g_coef(7,0)+b3*g_coef(8,0);
	r = 0+r*r_coef(0,0)+g*r_coef(1,0)+b*r_coef(2,0)+r2*r_coef(3,0)+g2*r_coef(4,0)+b2*r_coef(5,0)+r3*r_coef(6,0)+g3*r_coef(7,0)+b3*r_coef(8,0);

	//-- Combining channels and returning
	bgr[0] = b;
	bgr[1] = g;
	bgr[2] = r;
	Mat adjImage;
	merge(bgr,adjImage);
	return adjImage;
}

int main(int argc, char** argv){
	string img_type = string (argv[1]);
	bool bool_nir = img_type=="NIR";
	bool bool_vis = img_type=="VIS";
	bool bool_vis_CH = img_type=="VIS_CH";
	if(bool_vis | bool_vis_CH){
	//-- Putting black rectangle over the barcode area of input and background images
		Mat inputImage = imread(argv[2]);
		Mat adjBackground = imread(argv[3]);

	//-- Processing the VIS image
		Mat adjImage;
		float det=0;
		float D;
	//-- Color homography
		if(bool_vis_CH){
			MatrixXd rh, gh, bh;
			get_standardizations(inputImage, det, rh, gh, bh);
			adjImage = color_homography(inputImage,rh,gh,bh);
			D = 1-det;
		}else{
			adjImage = inputImage;
		}

	//-- Difference in images
		Mat dest;
		absdiff(adjBackground,adjImage,dest);
		vector<Mat> channels(3);
		split(dest,channels);
		Mat dest_blur;
		blur(channels[1], dest_blur, Size( 4, 4 ) );
		Mat dest_thresh;
		threshold(dest_blur,dest_thresh,40,255,THRESH_BINARY);
		Mat dest_dilate;
    	dilate(dest_thresh, dest_dilate, Mat(), Point(-1, -1), 5, 1, 1);
    	Mat dest_erode;
    	erode(dest_dilate,dest_erode, Mat(), Point(-1, -1), 5, 1, 1);

   	//-- Removing barcode
    	Mat lab;
    	cvtColor(adjImage, lab, CV_BGR2Lab);
    	vector<Mat> split_lab;
    	split(lab, split_lab);
    	Mat b_thresh1;
    	inRange(split_lab[2],90,141,b_thresh1);
    	Mat invSrc =  cv::Scalar::all(255) - b_thresh1;
    	Mat mask1;
    	bitwise_and(dest_erode,invSrc,mask1);

    //-- Remove edges of pot
    	Mat dest_lab;
    	cvtColor(dest, dest_lab, CV_BGR2Lab);
    	vector<Mat> channels_lab;
    	split(dest_lab, channels_lab);
    	Mat pot_thresh1;
    	inRange(channels_lab[2],0,120,pot_thresh1);
    	Mat pot_thresh2;
    	inRange(channels_lab[2],135,200,pot_thresh2);
    	Mat pot_or;
    	bitwise_or(pot_thresh1,pot_thresh2,pot_or);
    	Mat pot_dilate;
    	dilate(pot_or, pot_dilate, Mat(), Point(-1, -1), 2, 1, 1);
    	Mat pot_erode;
    	erode(pot_dilate,pot_erode, Mat(), Point(-1, -1), 3, 1, 1);
    	Mat pot_and;
    	bitwise_and(pot_erode,mask1,pot_and);
    	Mat pot_roi;
    	vector<Point> cc_pot = keep_roi(pot_and,Point(300,600),Point(1610,1310),pot_roi);

    //-- Remove blue stakes
    	Mat b_thresh;
    	inRange(split_lab[2],80,115,b_thresh);
    	Mat b_er;
    	erode(b_thresh,b_er, Mat(), Point(-1, -1), 1, 1, 1);
    	Mat b_roi;
    	vector<Point> cc1 = keep_roi(b_er,Point(300,600),Point(1610,1310),b_roi);
    	Mat b_dil;
    	dilate(b_roi,b_dil,Mat(),Point(-1, -1), 6, 1, 1);
    	Mat b_xor = pot_roi - b_dil;

    //-- ROI selector
    	Mat mask;
    	vector<Point> cc = keep_roi(b_xor,Point(550,0),Point(1810,1410),mask);

    //-- Getting numerical data
    	vector<double> shapes_data = get_shapes(cc,mask);
    	Mat hue_data = get_color(adjImage, mask);

    //-- Write shapes to file
    	string name_shape= string(argv[4]);
    	ofstream shape_file;
    	shape_file.open(name_shape.c_str(),ios_base::app);
    	shape_file << argv[2] << " ";
    	for(int i=0;i<20;i++){
    		shape_file << shapes_data[i];
    		if(i != 19){
    			shape_file << " ";
    		}
    	}
    	if(bool_vis_CH){
    		shape_file << " " << D;
    	}
    	shape_file << endl;
    	shape_file.close();

    //-- Write color to file
    	string name_hue= string(argv[5]);
    	ofstream hue_file;
    	hue_file.open(name_hue.c_str(),ios_base::app);
    	hue_file << argv[2] << " ";
    	for(int i=0;i<180;i++){
    		hue_file << hue_data.at<float>(i,0) << " ";
    	}
    	hue_file << endl;
    	hue_file.close();
	}
    //-- Processing the NIR image
	else if(bool_nir){
		//-- Read in image and background
    	Mat nirImage = imread(argv[2],0);
    	Mat nir_fixed = 1.591*nirImage-31.803;
    	namedWindow("Original",WINDOW_NORMAL);
    	    	        	    resizeWindow("Original",800,800);
    	    	        	    imshow("Original", nirImage);
    	namedWindow("Enhanced",WINDOW_NORMAL);
    	        	    resizeWindow("Enhanced",800,800);
    	        	    imshow("Enhanced", nir_fixed);
    	waitKey(0);
    	Mat nirBackground = imread(argv[3],0);

    	//-- Difference between image and background
		Mat dest_nir;
		absdiff(nirBackground,nirImage,dest_nir);
		Mat dest_nir_thresh;
		inRange(dest_nir,10,255,dest_nir_thresh);

		//-- Remove white stake
		Mat dest_stake;
		inRange(dest_nir,60,255,dest_stake);
		Mat dest_stake_dil;
		dilate(dest_stake, dest_stake_dil, Mat(), Point(-1, -1), 2, 1, 1);
		Mat kept_stake;
    	vector<Point> cc = keep_roi(dest_stake_dil,Point(270,183),Point(350,375),kept_stake);
    	Mat dest_sub = dest_nir_thresh - kept_stake;

        //-- ROI selector
    	Mat kept_mask_nir;
    	cc = keep_roi(dest_sub,Point(171,102),Point(470,363),kept_mask_nir);

        //-- Getting numerical data
    	Mat nir_data = get_nir(nirImage, kept_mask_nir);

        //-- Writing numerical data
    	string name_nir= string(argv[4]);
   		ofstream nir_file;
   		nir_file.open(name_nir.c_str(),ios_base::app);
   		nir_file << argv[2] << " ";
   		for(int i=0;i<255;i++){
   		   	nir_file << nir_data.at<float>(i,0) << " ";
   		}
   		nir_file << endl;
   		nir_file.close();
    }else{
    	cout << "First argument must be either VIS, VIS_CH, or NIR" << endl;
    }

	return 0;
}

/*
namedWindow("Image",WINDOW_NORMAL);
        	    resizeWindow("Image",800,800);
        	    imshow("Image", b_blur);
waitKey(0);
*/

