/**
@file utils.h
Application auxiliary functions

@brief Feature tracking using Lucas-Kanade algorithm

@date 18 May 2016

@author Wesley Serrano
*/

#include "CImg.h"
#include <ANN/ANN.h>          
#include <ANN/ANNx.h>         
#include <ANN/ANNperf.h>        
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <cassert>

using namespace cimg_library;
using namespace std;

typedef struct _point
{
  double x;
  double y;
} Point;

typedef struct _vector2
{
  double x;
  double y;
} Vector2;

typedef struct _vector5
{
  double x;
  double y;
  double L;
  double u;
  double v;
} Vector5;

double* convertRGBToLUV(double r, double g, double b);

typedef struct _color
{
  double r;
  double g;
  double b;

} Color;

Vector5 makeVec5(double x, double y, double r, double g, double b)
{
  Vector5 vec;

  vec.x = x;
  vec.y = y;

  double* Luv = convertRGBToLUV(r, g, b);

  vec.L = Luv[0];
  vec.u = Luv[1];
  vec.v = Luv[2];

  return vec;
}

Vector5 getColor(Vector5 v)
{
  Vector5 v2 = v;
  v2.x = 0.0;
  v2.y = 0.0;

  return v2;
}

Vector5 getSpace(Vector5 v)
{
  Vector5 v2 = v;
  v2.L = 0.0;
  v2.u = 0.0;
  v2.v = 0.0;

  return v2;
}

Vector5 add(Vector5 v1, Vector5 v2)
{
  Vector5 v;

  v.x = v1.x + v2.x;
  v.y = v1.y + v2.y;
  v.L = v1.L + v2.L;
  v.u = v1.u + v2.u;
  v.v = v1.v + v2.v;

  return v;
}

Vector5 divide(Vector5 v, const double c)
{
  Vector5 v1;

  if(c == 0.0) return makeVec5(0,0,0,0,0);
  v1.x = v.x/c;
  v1.y = v.y/c;
  v1.L = v.L/c;
  v1.u = v.u/c;
  v1.v = v.v/c;

  return v1;
}

Vector5 multiply(Vector5 v, const double c)
{
  Vector5 v1;

  v1.x = c*v.x;
  v1.y = c*v.y;
  v1.L = c*v.L;
  v1.u = c*v.u;
  v1.v = c*v.v;

  return v1;
}

Vector5 subtract(Vector5 v1, Vector5 v2)
{
  Vector5 v;

  v.x = v1.x - v2.x;
  v.y = v1.y - v2.y;
  v.L = v1.L - v2.L;
  v.u = v1.u - v2.u;
  v.v = v1.v - v2.v;

  return v;
}

double length(Vector5 v)
{
  return sqrt(v.x*v.x + v.y*v.y + v.L*v.L + v.u*v.u + v.v*v.v);
}

double* convertRGBToLUV(double r, double g, double b)
{
  r /= 255.0; g /= 255.0; b /= 255.0;
  const double X = 0.412453*r + 0.35758*g + 0.180423*b, Xr = 0.95047;
  const double Y = 0.212671*r + 0.71516*g + 0.072169*b, Yr = 1.0;
  const double Z = 0.019334*r + 0.119193*g + 0.950227*b, Zr = 1.08883;

  const double y = (Y/Yr), epsilon = 0.008856;

  const double L = (y > epsilon)? (116*pow(y,1/3) - 16):903.3*y;

  const double denominator = X + 15*Y + 3*Z;

  const double u0 = (denominator != 0) ? (4*X)/denominator : 4, ur = 0.19784977571475;
  const double v0 = (denominator != 0) ? (9*Y)/denominator : 9.0/15.0, vr = 0.46834507665248;

  const double u = 13.0*L*(u0 - ur);
  const double v = 13.0*L*(v0 - vr);

  double* Luv = new double[3];

  Luv[0] = L;
  Luv[1] = u;
  Luv[2] = v;

  return Luv;
}

/**
  Returns the maximum of three values of the colors channels of a pixel

  @param r The value of the red channel
  @param g The value of the green channel
  @param b The value of the blue channel

  @return The maximum of the color channels
*/
double MaxTone(double r,double g,double b)
{
  double aux = max(r,g);
  return max(aux,b);
}

/**
  Returns the minimum of three values of the colors channels of a pixel

  @param r The value of the red channel
  @param g The value of the green channel
  @param b The value of the blue channel

  @return The minimum of the color channels
*/
float MinTone(float r,float g,float b)
{
  float aux = min(r,g);
  return min(aux,b);
}

/**
  Returns the average of three values of the colors channels of a pixel

  @param r The value of the red channel
  @param g The value of the green channel
  @param b The value of the blue channel

  @return The average of the color channels
*/
float AvgTone(float r,float g,float b)
{
  return (r + g + b)/3.0;
}

/**
  Creates a window to display a image given as input parameter

  @param image Image to be shown on a window
*/

void displayImage(CImg<double> image)
{
  CImgDisplay display(image, "Image");

  while (!display.is_closed()) 
  {
    display.wait();
  }
}

/**
  Creates a window to display a image given as input parameter and draws a rectangle on the giver coordinates

  @param image Image to be shown on a window
*/

void displayImageWithRectangle(CImg<double> image, int x0, int y0, int x1, int y1)
{
  double color[] = {0, 255, 0};
  CImgDisplay display(image.draw_line(x0, y0, x0, y1, color)
                           .draw_line(x0, y1, x1, y1, color)
                           .draw_line(x1, y1, x1, y0, color)
                           .draw_line(x0, y0, x1, y0, color), "Image");

  while (!display.is_closed()) 
  {
    display.wait();
  }
}

vector<Point> getBoundary(CImg<double> image)
{
  double color[] = {0, 255, 0};
  int x0 = 0, y0 = 0, x1 = image.width() - 1, y1 = image.height() - 1;
  CImg<double> displayImage = image;
  CImgDisplay display(displayImage.draw_line(x0, y0, x0, y1, color)
                           .draw_line(x0, y1, x1, y1, color)
                           .draw_line(x1, y1, x1, y0, color)
                           .draw_line(x0, y0, x1, y0, color), "Image");

  while (!display.is_closed()) 
  {
    if(display.button()&1)
    {
      x0 = display.mouse_x();
      y0 = display.mouse_y();
      
      displayImage = image;
      displayImage.draw_line(x0, y0, x0, y1, color)
                  .draw_line(x0, y1, x1, y1, color)
                  .draw_line(x1, y1, x1, y0, color)
                  .draw_line(x0, y0, x1, y0, color).display(display);
    }
    if(display.button()&2)
    {
      x1 = display.mouse_x();
      y1 = display.mouse_y();
      displayImage = image;
      displayImage.draw_line(x0, y0, x0, y1, color)
                  .draw_line(x0, y1, x1, y1, color)
                  .draw_line(x1, y1, x1, y0, color)
                  .draw_line(x0, y0, x1, y0, color).display(display);
    }
    display.wait();
  }

  Point firstBoundaryPoint, secondBoundaryPoint;
  firstBoundaryPoint.x = x0;
  firstBoundaryPoint.y = y0;
  secondBoundaryPoint.x = x1;
  secondBoundaryPoint.y = y1;
  
  vector<Point> boundary;
  boundary.push_back(firstBoundaryPoint);
  boundary.push_back(secondBoundaryPoint);

  return boundary;
}

/**
  Creates a window to display a vector of images given as input parameter.
  Displays one image at a time

  @param images The vector containing the images to be shown on a window
*/
void displayImages(vector< CImg<double> > images)
{
  for(int i = 0; i < images.size(); i++) displayImage(images[i]);
}

/**
  Creates a gray scale version of an image

  @param image Image to be transformed
  @param scale Indicates whether or not the RGB channels values should be in the range [0, 1] or [0, 255] (as integer values)

  @return Image in gray scale
*/
CImg<double> makeImageGray(CImg<double> image, bool scale = true)
{
  const bool scaleFactor = scale ? 255 : 1;
  CImg<double> result = image;

  cimg_forXY(result, x, y)
  {
    if(image.spectrum() > 1)
    {
      const double r = result(x, y, 0, 0)/scaleFactor;
      const double g = result(x, y, 0, 1)/scaleFactor;
      const double b = result(x, y, 0, 2)/scaleFactor;

      result(x, y, 0, 0) = 0.299 * r + 0.587 * g + 0.114 * b;
      result(x, y, 0, 1) = 0.299 * r + 0.587 * g + 0.114 * b;
      result(x, y, 0, 2) = 0.299 * r + 0.587 * g + 0.114 * b;
    }
    else
    {
      const double c = result(x, y, 0, 0)/scaleFactor;
      result(x, y, 0, 0) = c;
    }
  }

  return result;
}


/**
  Creates a window to display a image given as input parameter
  with green dots on the points given as the other parameter

  @param image Image to be marked and shown on a window
  @param pointsToMark The coordinates of the points to be marked on the image
*/
void markPointsOnImage(CImg<double> image, vector<Point> pointsToMark)
{
  const int WIDTH = image.width(), HEIGHT = image.height();
  CImg<double> points(WIDTH,HEIGHT,1,3,0);

   cimg_forXY(image,x,y)
   {
     const double c = image(x, y, 0, 0);

    points(x, y, 0, 0) = c;
    points(x, y, 0, 1) = c;
    points(x, y, 0, 2) = c;
   }

  for(vector<Point>::iterator it = pointsToMark.begin(); it != pointsToMark.end(); ++it)
  {
    const double x = it->x;
    const double y = it->y;

    if(x < 0 || x >= WIDTH) continue;
    if(y < 0 || y >= HEIGHT) continue;

    points(x, y, 0, 0) = 0;
    points(x, y, 0, 1) = 255;
    points(x, y, 0, 2) = 0;
  }

  displayImage(points);
}

/**
  Sums the rgb channels of two images and save the result on another image.
  If one of the values is greater than 255, it is changed to 255.

  @param image1 First image to be used
  @param image2 Second image to be used 

  @return Image containing the result of the sum
*/
CImg<double> add(CImg<double> image1, CImg<double> image2)
{
   const int WIDTH1 = image1.width(), HEIGHT1 = image1.height();
   const int WIDTH2 = image2.width(), HEIGHT2 = image2.height();

   //assert(WIDTH1 == WIDTH2 && HEIGHT1 == HEIGHT2);

   CImg<double> newImage(WIDTH1,HEIGHT1,1,3,0);
   
   cimg_forXY(image1,x,y)
   {
     const double r1 = image1(x,y,0,0), g1 = image1(x,y,0,1), b1 = image1(x,y,0,2);
     const double r2 = image2(x,y,0,0), g2 = image2(x,y,0,1), b2 = image2(x,y,0,2);

     newImage(x,y,0,0) = r1 + r2;
     newImage(x,y,0,1) = g1 + g2;
     newImage(x,y,0,2) = b1 + b2;
   }

   return newImage;
}

/**
  Subtracts the rgb channels of two images and save the result on another image.
  For every color channel of every pixel, the absolute value of the difference is saved. 

  @param image1 First image to be used
  @param image2 Second image to be used 

  @return Image containing the result of the subtraction
*/
CImg<double> subtract(CImg<double> image1, CImg<double> image2)
{
   const int WIDTH1 = image1.width(), HEIGHT1 = image1.height();
   const int WIDTH2 = image2.width(), HEIGHT2 = image2.height();

   //assert(WIDTH1 == WIDTH2 && HEIGHT1 == HEIGHT2);


   CImg<double> newImage(WIDTH1,HEIGHT1,1,3,0);
   
   cimg_forXY(image1,x,y)
   {
     const double r1 = image1(x,y,0,0), g1 = image1(x,y,0,1), b1 = image1(x,y,0,2);
     const double r2 = image2(x,y,0,0), g2 = image2(x,y,0,1), b2 = image2(x,y,0,2);

     newImage(x,y,0,0) = abs(r1 - r2);
     newImage(x,y,0,1) = abs(g1 - g2);
     newImage(x,y,0,2) = abs(b1 - b2);
   }

   return newImage;
}

CImg<double> unidimensionalConvolution(CImg<double> image, double* convolutionMatrix, const int matrixSize, const bool vertical, const bool grayScale = false)
{
   const int offset = (int) matrixSize/2;
   const int WIDTH = image.width(), HEIGHT = image.height();
   CImg<double> resultImage(WIDTH,HEIGHT,1, grayScale ? 1 : 3,0);

   cimg_forXY(image,x,y)
   {
     double imageR, imageG, imageB;
     double imageC;

     double r = 0.0, g = 0.0, b = 0.0;
     double c = 0.0;

     for(int i = -offset; i <= offset; i++)
     {
      const int pixelY = vertical ? y + i : y;
      const int matrixCell = i + offset;
      bool pixelOutOfImage;
      
      const int pixelX = vertical? x : x + i;

      if((pixelY < 0 || pixelY >= HEIGHT) || (pixelX < 0 || pixelX >= WIDTH)) pixelOutOfImage = true;
      else pixelOutOfImage = false;

      if(grayScale)
      {
        imageC = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,0); 
        c += imageC * convolutionMatrix[matrixCell];
      }
      else
      {
       imageR = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,0); 
       imageG = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,1); 
       imageB = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,2); 

       r += imageR * convolutionMatrix[matrixCell];
       g += imageG * convolutionMatrix[matrixCell];
       b += imageB * convolutionMatrix[matrixCell];
      }
      
     }

     if(grayScale) resultImage(x,y,0,0) = c;
     else
     { 
      resultImage(x,y,0,0) = r;
      resultImage(x,y,0,1) = g;
      resultImage(x,y,0,2) = b;
     }
   }   

   return resultImage;
}

/**
  Convolves an image with a matrix
  @param image Image to apply to convolution
  @param convolutionMatrix The convolution matrix to be applied. Must be square.
  @param matrixSize The square matrix dimension
  @param grayScale Indicates Whether or not the image to be convolved is a gray scale image

  @return The convolved image
*/
CImg<double> convolve(CImg<double> image, double** convolutionMatrix, const int matrixSize, const bool grayScale = false)
{
   const int offset = (int) matrixSize/2;
   const int WIDTH = image.width(), HEIGHT = image.height();
   CImg<double> resultImage(WIDTH,HEIGHT,1, grayScale ? 1 : 3,0);

   cimg_forXY(image,x,y)
   {
     double imageR, imageG, imageB;
     double imageC;

     double r = 0.0, g = 0.0, b = 0.0;
     double c = 0.0;

     for(int i = -offset; i <= offset; i++)
     {
      const int pixelY = y + i;
      const int row = i + offset;
      bool pixelOutOfImage;

      for(int j = -offset; j <= offset; j++)
      {
        const int pixelX = x + j;
        const int column = j + offset;

        if((pixelY < 0 || pixelY >= HEIGHT) || (pixelX < 0 || pixelX >= WIDTH)) pixelOutOfImage = true;
        else pixelOutOfImage = false;

        if(grayScale)
        {
          imageC = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,0); 
          c += imageC * convolutionMatrix[row][column];
        }
        else
        {
         imageR = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,0); 
         imageG = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,1); 
         imageB = pixelOutOfImage ? 0.0 : image(pixelX,pixelY,0,2); 

         r += imageR * convolutionMatrix[row][column];
         g += imageG * convolutionMatrix[row][column];
         b += imageB * convolutionMatrix[row][column];
        }
      }
     }

     if(grayScale) resultImage(x,y,0,0) = c;
     else
     { 
      resultImage(x,y,0,0) = r;
      resultImage(x,y,0,1) = g;
      resultImage(x,y,0,2) = b;
     }
   }   

   return resultImage;
}

/**
  Creates a Gaussian filter matrix
  @param filterSize The matrix dimension. Must be odd
  @param delta The delta parameter of the gaussian function
  @param correction Indicates whether or not the matrix values should be adjusted

  @return The filter matrix
*/
double** makeGaussianFilter(const int filterSize, const double delta, const bool correction = false)
{
   assert(filterSize%2 == 1);

   double** filterMatrix;
   filterMatrix = new double*[filterSize];
   for (int i = 0; i < filterSize; i++) filterMatrix[i] = new double[filterSize];

   double sum = 0.0;

  for (int i = 0; i < filterSize; i++)
  {
    for (int j = 0; j < filterSize; j++)
    {
      int i2 = i-(filterSize-2);
      int j2 = j-(filterSize-2);
      double value = (exp((-(i2*i2+j2*j2))/(2*delta*delta)))/(2*M_PI*delta*delta);
      filterMatrix[j][i] = value;
      sum += value;
    }
  }

  if(correction) sum = sum/2.0;

  for (int i = 0; i < filterSize; i++)
  {
    for (int j = 0; j < filterSize; j++)
    {
      double value2 = filterMatrix[j][i]/sum;
      filterMatrix[j][i] = value2;
    }
  }

   return filterMatrix;
}

double** makeConstantGaussianFilter(const double correction)
{
   double** filterMatrix;
   filterMatrix = new double*[5];
   for (int i = 0; i < 5; i++) filterMatrix[i] = new double[5];
   double filterWeight = correction ? 8.0 : 16.0;
    filterWeight *= filterWeight;

   filterMatrix[0][0] = 1.0/filterWeight; filterMatrix[0][1] = 4.0/filterWeight; filterMatrix[0][2] = 6.0/filterWeight; filterMatrix[0][3] = 4.0/filterWeight; filterMatrix[0][4] = 1.0/filterWeight;
   filterMatrix[1][0] = 4.0/filterWeight; filterMatrix[1][1] = 16.0/filterWeight; filterMatrix[1][2] = 24.0/filterWeight; filterMatrix[1][3] = 16.0/filterWeight; filterMatrix[1][4] = 4.0/filterWeight;
   filterMatrix[2][0] = 6.0/filterWeight; filterMatrix[2][1] = 24.0/filterWeight; filterMatrix[2][2] = 36.0/filterWeight; filterMatrix[2][3] = 24.0/filterWeight; filterMatrix[2][4] = 6.0/filterWeight;
   filterMatrix[3][0] = 4.0/filterWeight; filterMatrix[3][1] = 16.0/filterWeight; filterMatrix[3][2] = 24.0/filterWeight; filterMatrix[3][3] = 16.0/filterWeight; filterMatrix[3][4] = 4.0/filterWeight;
   filterMatrix[2][0] = 1.0/filterWeight; filterMatrix[4][1] = 4.0/filterWeight; filterMatrix[4][2] = 6.0/filterWeight; filterMatrix[4][3] = 4.0/filterWeight; filterMatrix[4][4] = 1.0/filterWeight;

   return filterMatrix;
}

double* makeConstantUnidimensionalGaussianFilter(const double correction)
{
   double* filterMatrix;
   filterMatrix = new double[5];

   double filterWeight = correction ? 8.0 : 16.0;

   filterMatrix[0] = 1.0/filterWeight; filterMatrix[1] = 4.0/filterWeight; filterMatrix[2] = 6.0/filterWeight; filterMatrix[3] = 4.0/filterWeight; filterMatrix[4] = 1.0/filterWeight;

   return filterMatrix;
}

/**
  Applies a 3x3 Gaussian filter on a image
  @param image image to apply the filter
  @param correction Wether a correction factor should or not be applied to the result
  @param delta The delta parameter in the gaussian filter generation

  @return The original image blurred
*/

CImg<double> blur3x3(CImg<double> image, const bool correction = false, const double delta = 1.5, const bool grayScale = false)
{
   const int filterSize = 3;

   double** convolutionMatrix;
   convolutionMatrix = makeGaussianFilter(filterSize, delta);
   
   return convolve(image, convolutionMatrix, filterSize, grayScale);
}

/**
  Applies a 5x5 Gaussian filter on a image
  @param image image to apply the filter
  @param correction Wether a correction factor should or not be applied to the result
  @param delta The delta parameter in the gaussian filter generation

  @return The original image blurred
*/

CImg<double> blur5x5(CImg<double> image, const bool correction = false, const double delta = 1.5, const bool grayScale = false)
{   
   const int filterSize = 5;

   double** convolutionMatrix;
   convolutionMatrix = makeGaussianFilter(filterSize, delta, correction);

   //convolutionMatrix = makeConstantGaussianFilter(correction);

   return convolve(image, convolutionMatrix, filterSize, grayScale);
   
   /*double* convolutionMatrix;
   convolutionMatrix = makeConstantUnidimensionalGaussianFilter(correction);

   CImg<double> verticallyFiltered = unidimensionalConvolution(image, convolutionMatrix, filterSize, true, grayScale);

   CImg<double> horizontallyFiltered = unidimensionalConvolution(verticallyFiltered, convolutionMatrix, filterSize, false, grayScale);
   

   return horizontallyFiltered;*/
}

/**
  Applies a Gaussian filter on a image. The filter matrix size depends on a parameter. If the filter size is not 3 or 5, it does nothing.
  @param image The image to apply the filter
  @param filterSize The filter matrix size.
  @param correction Wether a correction factor should or not be applied to the result
  @param delta The delta parameter in the gaussian filter generation

  @return The original image blurred
*/
CImg<double> blur(CImg<double> image, const int filterSize, const bool correction = false, const double delta = 2.0, const bool grayScale = false)
{
  if(filterSize == 3) return blur3x3(image, correction, delta, grayScale);
  else if(filterSize == 5) return blur5x5(image, correction, delta, grayScale);
}

/**
  Double
  @param image First image to be expanded

  @return The expanded image
*/
CImg<double> expandImage(CImg<double> image)
{
   const int WIDTH = image.width(), HEIGHT = image.height();
   
   CImg<double> expandedImage(2*WIDTH, 2*HEIGHT,1,3,0);

   cimg_forXY(expandedImage,x,y)
   {
     if(x%2 == 0 && y%2 == 0)
     {
       const int px = x/2, py = y/2;
       double r = image(px,py,0,0), g = image(px,py,0,1), b = image(px,py,0,2);   
       
       expandedImage(x,y,0,0) = r;
       expandedImage(x,y,0,1) = g;
       expandedImage(x,y,0,2) = b;
     }
   }

   return expandedImage;
}

/**
  Reduces a image by half
  @param image First image to be reduced

  @return The reduced image
*/
CImg<double> reduceImage(CImg<double> image, const bool grayScale = false)
{
   const int WIDTH = image.width(), HEIGHT = image.height();
   
   CImg<double> reducedImage(WIDTH%2 == 0 ? WIDTH/2 : (WIDTH + 1)/2,HEIGHT%2 == 0 ? HEIGHT/2 : (HEIGHT + 1)/2,1, grayScale ? 1 : 3,0);

   cimg_forXY(image,x,y)
   {
     if(x%2 == 0 && y%2 == 0)
     {
       if(grayScale)
       {
        double c = image(x,y,0,0);
        reducedImage(x/2,y/2,0,0) = c;        
       }
       else
       {
        double r = image(x,y,0,0), g = image(x,y,0,1), b = image(x,y,0,2);   
       
        reducedImage(x/2,y/2,0,0) = r;
        reducedImage(x/2,y/2,0,1) = g;
        reducedImage(x/2,y/2,0,2) = b;
       }
     }
   }

   return reducedImage;
}

/**
  Makes a gaussian pyramid of images from an original image. 
  
  @param imageFilePath The path to the image
  @param numberOfImages Output parameter to save the number of generated images from the original
  @param gaussianFilterSize The size of the gaussian filter to be used

  @return Image containing all the pyramid images 
*/
vector< CImg<double> > makeGaussianPyramid(CImg<double> image, const int gaussianFilterSize, const bool grayScale = false)
{
  vector< CImg<double> > pyramid;

  const int MINIMUM_WIDTH = 32;
  const double DELTA = 2.0;
  const int width = image.width(), height = image.height();
    
  const int numberOfImages = (int)(log2(width) - log2(MINIMUM_WIDTH));
  
  CImg<double> firstBlur = blur(image, gaussianFilterSize, false, DELTA, grayScale);
  pyramid.push_back(grayScale ? firstBlur : image);

  //CImg<double> blurredImage = blur(pyramid[0], gaussianFilterSize, false, DELTA, grayScale);

  //CImg<double> silhouetteImage = subtract(pyramid[0], blurredImage);

  for(int i = 1; i < numberOfImages; i++)
  { 
    CImg<double>  blurredImage = blur(pyramid[i - 1], gaussianFilterSize, false, DELTA, grayScale);
    
    pyramid.push_back(reduceImage(blurredImage, grayScale));
  }

  return pyramid;
}
