#include <iostream>
#include <opencv2/opencv.hpp>
#include "include/ImageTransformation.h"
#include <time.h>
#include <mpi.h>

using namespace cv;
using namespace std;

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cerr << "Invalid arguments number" << endl;
        return -1;
    }

    Mat originalImage = imread(argv[1]);

    int imageHeight = originalImage.rows;
    int imageWidth = originalImage.cols;

//    Mat medianedImage(imageHeight, imageWidth, CV_8UC3);
//    Mat averagedImage(imageHeight, imageWidth, CV_8UC3);
    Mat bilateralImage(imageHeight, imageWidth, CV_8UC3);

    if (originalImage.empty() || originalImage.depth() != CV_8U)
    {
        cerr << "Image " << argv[0] << " is empty or type invalid" << endl;
        return -1;
    }

    Mat bgr[3];
    split(originalImage, bgr);


//OpenCV uses BGR color order!!!
    auto *B = new ImageDataClass(bgr[0].data, imageHeight, imageWidth, 1);
    auto *G = new ImageDataClass(bgr[1].data, imageHeight, imageWidth, 1);
    auto *R = new ImageDataClass(bgr[2].data, imageHeight, imageWidth, 1);

    double spatialSigma = 24;
    double rSigma = 1;
    int kernelSize = 9;

    bgr[0].datastart = bilateral_filer(B, kernelSize, spatialSigma, rSigma, NN_CLONE);
    bgr[1].datastart = bilateral_filer(G, kernelSize, spatialSigma, rSigma, NN_CLONE);
    bgr[2].datastart = bilateral_filer(R, kernelSize, spatialSigma, rSigma, NN_CLONE);

    merge(bgr, 3, bilateralImage);

//    delete B;
//    delete G;
//    delete R;
//
//    bgr[0].data = convolve<unsigned char>(B, kernelSize, median_filter, nullptr, NN_CLONE);
//    bgr[1].data = convolve<unsigned char>(G, kernelSize, median_filter, nullptr, NN_CLONE);
//    bgr[2].data = convolve<unsigned char>(R, kernelSize, median_filter, nullptr, NN_CLONE);
//    merge(bgr, 3, medianedImage);
//
//    Mat average_bgr[3];
//    split(originalImage, average_bgr);
//
//OpenCV uses BGR color order!!!
//    ImageDataClass aB(average_bgr[0].data, imageHeight, imageWidth, 1);
//    ImageDataClass aG(average_bgr[1].data, imageHeight, imageWidth, 1);
//    ImageDataClass aR(average_bgr[2].data, imageHeight, imageWidth, 1);
//
//    int kernelSize = 11;
//    auto *hat_kernel_matrix = new TransformationMatrix<double>(kernelSize, multiply_each, hat_kernel);
//
//    average_bgr[0].data = convolve(&aB, kernelSize, sum_pixel_values, hat_kernel_matrix, NN_CLONE);
//    average_bgr[1].data = convolve(&aG, kernelSize, sum_pixel_values, hat_kernel_matrix, NN_CLONE);
//    average_bgr[2].data = convolve(&aR, kernelSize, sum_pixel_values, hat_kernel_matrix, NN_CLONE);
//
//    delete hat_kernel_matrix;
//
//    merge(average_bgr, 3, averagedImage);

//    imwrite("leoMedian11.jpg", medianedImage);
//    imwrite("leoAvg11.jpg", averagedImage);
    imwrite("leoBilat_r1_s_24_9.jpg", bilateralImage);
    namedWindow("original");
//    namedWindow("medianFilter");
//    namedWindow("averageHat");
    namedWindow("bilateral");
    imshow("original", originalImage);
//    imshow("medianFilter", medianedImage);
//    imshow("averageHat", averagedImage);
    imshow("bilateral", bilateralImage);
//    imwrite(argv[2], bilateralImage);
    waitKey(0);
    destroyAllWindows();


}


