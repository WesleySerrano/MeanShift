#include "utils.h"

vector<Vector5> initializeKernels(CImg<double> image)
{
  vector<Vector5> result;
  cimg_forXY(image, x, y)
  {
    Vector5 vec = makeVec5(x, y, image(x, y, 0, 0), image(x, y, 0, 1), image(x, y, 0, 2));

    result.push_back(vec);
  }

  return result;
}

double makeKernelProfile(Vector5 v, const int h)
{
  const double vectorLength = length(v);

  return exp(-0.5*(vectorLength * vectorLength)/(h*h));
}

double makeKernelConstant(const int hr, const int hs, const int d, const int n)
{
  const int p = d - 2;
  const double pow1 = pow(2*M_PI, d/2.0), pow2 = pow(hr, p), pow3 = pow(hs, 2.0);

  return 1.0/(pow1 * pow2 * pow3 * n);
}

double meanShiftExponential(Vector5 v, const int hr, const int hs)
{
  Vector5 colorVector = getColor(v), spaceVector = getSpace(v);
  Vector5 scaledColor = divide(colorVector, hr), scaledSpace = divide(spaceVector, hs);
  const double colorLength = length(scaledColor);
  const double spaceLength = length(scaledSpace);

  const double colorWeight = exp(-0.5 * (colorLength * colorLength));
  const double spaceWeight = exp(-0.5 * (spaceLength * spaceLength));

  return colorWeight * spaceWeight;
}

ANNkd_tree* makeKdTree(vector<Vector5> kernels, const int d)
{
  const int SIZE = kernels.size();

  ANNpointArray points = annAllocPts(SIZE, d);
 
  for(int i = 0; i < SIZE; i++)
  {
    Vector5 kernel = kernels[i];
    ANNpoint point = annAllocPt(d);

     point[0] = kernel.x;
     point[1] = kernel.y;
     point[2] = kernel.L;
     point[3] = kernel.u;
     point[4] = kernel.v;

     points[i] = point;
  }

  return new ANNkd_tree(points, SIZE, d);
}

vector<Vector5> meanShiftFiltering(vector<Vector5> kernels, const int hr, const int hs)
{
 vector<Vector5> filteredKernels = kernels;
 vector<bool> stopped;
 const int d = 5, MAX_NEIGHBORS = 30;
 const double epsilon = 0.00000001;

 ANNkd_tree* kdTree;
 kdTree = makeKdTree(kernels, d);

 for(int i = 0; i < filteredKernels.size(); i++)
 {
   ANNpoint kernelPoint = annAllocPt(d);
   Vector5 kernel = kernels[i];
   Vector5 Yt = makeVec5(0,0,0,0,0);

   kernelPoint[0] = kernel.x;
   kernelPoint[1] = kernel.y;
   kernelPoint[2] = kernel.L;
   kernelPoint[3] = kernel.u;
   kernelPoint[4] = kernel.v;

   ANNidxArray neighborhood; neighborhood = new ANNidx[MAX_NEIGHBORS];
   ANNdistArray distances; distances = new ANNdist[MAX_NEIGHBORS];
   double distance = 1.0;
   while(distance > epsilon)
   {
     Vector5 numerator = makeVec5(0,0,0,0,0);
     double denominator = 0.0;

     kdTree->annkSearch(kernelPoint, MAX_NEIGHBORS, neighborhood, distances);

     for(int j = 0; j < MAX_NEIGHBORS; j++)
     {
        Vector5 neighbor = kernels[neighborhood[j]];
        Vector5 difference = subtract(neighbor, kernel);

        const double exponentialWeight = meanShiftExponential(difference, hr, hs);
        denominator += exponentialWeight;
        numerator = add(numerator, multiply(neighbor, exponentialWeight));
     }

     Yt = divide(numerator, denominator);
     Vector5 displacementVector = subtract(Yt, kernel);

     kernel = Yt;

     kernelPoint[0] = kernel.x;
     kernelPoint[1] = kernel.y;
     kernelPoint[2] = kernel.L;
     kernelPoint[3] = kernel.u;
     kernelPoint[4] = kernel.v;

     distance = length(displacementVector);
   }

   filteredKernels[i] = Yt;
 } 

 return filteredKernels;
}

CImg<double> makeFilteredImage(CImg<double> image, vector<Vector5> originalKernels, vector<Vector5> filteredKernels)
{
  CImg<double> result = image;

  for(int i = 0; i < filteredKernels.size(); i++)
  {
    Vector5 filteredKernel = filteredKernels[i];
    Vector5 originalKernel = originalKernels[i];

    const int x0 = originalKernel.x, y0 = originalKernel.y;
    const int x1 = filteredKernel.x, y1 = filteredKernel.y;

    result(x0, y0, 0, 0) = image(x1, y1, 0, 0);
    result(x0, y0, 0, 1) = image(x1, y1, 0, 1);
    result(x0, y0, 0, 2) = image(x1, y1, 0, 2);
  }

  return result;
}

int main (int argc, char **argv)
{
   if(argc < 4)
   {
   	 cout << "Insufficient arguments\n";
     cout << "You must type: " << argv[0] << " <image> <hs value> <hr value>";
   	 exit(-1);
   }

  CImg<double> image(argv[1]);
  const int hs = atoi(argv[2]), hr = atoi(argv[3]);

  displayImage(image);

  clock_t begin_time = clock();
  vector<Vector5> firstKernels = initializeKernels(image);
  cout << "Filtering\n";
  vector<Vector5> newPoints = meanShiftFiltering(firstKernels, hr, hs);
  cout << "Finished in " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds\n";
  CImg<double> filteredImage = makeFilteredImage(image, firstKernels, newPoints);
  displayImage(filteredImage);

  filteredImage.save("filtered.png");

  return 0;
}