#define cimg_display 0
#include "CImg.h"
#include <Eigen/Dense>
#include <cstring>
#include <vector>
#include <iostream>


using namespace cimg_library;
using namespace std;



void seamer(vector<Eigen::Vector3d> &image, int &height, int &width){
	
  double MaxEnergy = sqrtf(7500.0); //50^2 *3

  double energy;
  //int width= input.width();
  //int height= input.height();
  vector<double> energies;
  vector<double> costs;

  vector<Eigen::Vector2d> gradient;
  gradient.resize(width * height);
  
  energies.resize(gradient.size());
  
  vector<int> seam;
  seam.resize(height);
  
  costs.resize(energies.size());
  
  //compute gradient
  for(int i = 0; i < height; i++){
	  for(int j = 0; j < width; j++){
		//-----------x stuff------------
		//starting at the left (the first pixel)
		if(j == 0){
		//x component
		gradient[i * width + j][0] = (image[i* width + j + 1] - image[i* width + j]).norm()/2;
		}
		//at the right side of the image (the last pixel)
		else if(j == width - 1){
		//x component
		gradient[i * width + j][0] = (image[i* width + j] - image[i* width + j - 1]).norm()/2;
		}
		else{
		//x component
		gradient[i * width + j][0] = (image[i* width + j + 1] - image[i* width + j -1]).norm()/2;
		}
		
		//----------y stuff------------
		//starting at the top
		if(i == 0){
		//y component
		gradient[i * width + j][1] = (image[i* width + j + width] - image[i* width + j]).norm()/2;	
		}
		//at the bottom of the image
		else if(i == height - 1){
		//y component
		gradient[i * width + j][1] = (image[i* width + j] - image[i* width + j - width]).norm()/2;		  
		}
		//everywhere else
		else{
		//y component
		gradient[i * width + j][1] = (image[i* width + j + width] - image[i* width + j - width]).norm()/2;		  
		}
		
	  }
  }
  //energy function
  for(int i = 0; i < gradient.size(); i++){
	  energies[i] = sqrtf(gradient[i][0] * gradient[i][0] + gradient[i][1] * gradient[i][1]);
  }
  
  //compute cost of seam
  for(int i = 0; i < height; i++){
	  for(int j = 0; j < width; j++){
		  //first row
		  if(i == 0){
			costs[i * width + j] = energies[i * width + j];
			
		  }
		  //left edge
		  else if(j == 0){
			  costs[i * width + j] = energies[i * width + j] + min(costs[(i- 1) * width + (j)],costs[(i-1) * width + (j+1)]);
					 
		 }
		  //right edge
		  else if(j == width - 1){
			  costs[i * width + j] = energies[i * width + j] + min(costs[(i- 1) * width + (j-1)],costs[(i- 1) * width + (j)]);
			
		 }
		  
		  //second row to the last row
		  else{
			costs[i * width + j] = energies[i * width + j] + min(min(costs[(i- 1) * width + (j-1)],costs[(i- 1) * width + (j)]),costs[(i-1) * width + (j+1)]);
			
		}
		
	  }
  }
  
  //create the seam
  seam[height - 1] = (height-1) * width;
  
  for(int j = 0; j < width; j++){
	if(costs[(height-1) * width + j] < costs[seam[height - 1]]){
		
		seam[height-1] = (height-1) * width + j;
		}	 
	} 
	//column of seam
	seam[height-1] = seam[height-1]%width;
  
  
  //based on min cost
	for(int i = height - 2; i > -1; i--){
		seam[i] = seam[i + 1];//sets it to prev column
		//first check
		//if(i == height - 2){
		double minCost = 1920*1920;
		//int upperb, lowerb; //upper and lower bounds
		
		if(seam[i+1] == 0){
			if(seam[i+1] == width-1){
				for(int j = 0; j < width; j++){
					if(costs[i*width +j] < minCost){
						minCost = costs[i *width + j];
						seam[i] = j;
					}
				}
			}
			else{
				for(int j = 0; j < seam[i+1]+2; j++){
					if(costs[i*width +j] < minCost){
						minCost = costs[i *width + j];
						seam[i] = j;
					}
				}
			}
			
		}
		else{
			if(seam[i+1] == width-1){
				for(int j = seam[i+1]-1; j < width; j++){
					if(costs[i*width +j] < minCost){
						minCost = costs[i *width + j];
						seam[i] = j;
					}
				}
			}
			else{
				for(int j = seam[i+1]-1; j < seam[i+1]+2; j++){
					if(costs[i*width +j] < minCost){
						minCost = costs[i *width + j];
						seam[i] = j;
					}
				}
			}
			
		}
	}
		
 
 //-----resize the image----------
 //include every pixel that isnt in the seam
 vector<Eigen::Vector3d> preImage;
 //preImage.resize(width*height);
 
 preImage = image;
 
 int index = 0;
 int nRemoved = 0;
 
 
 image.resize((width-1)*height);
 
 for(int i = 0; i < height; i++){
	 for(int j = 0; j < width; j++){
		if(j != seam[i]){
			image[index] = preImage[i*width+j];
			index++;
		 }
	 }
 }
	
}




int main(int argc, char *argv[]) {
  CImg<double> input(argv[1]);
  CImg<double> lab = input;
  lab.RGBtoLab();
  
  vector<Eigen::Vector3d> image;
  image.resize(input.width() * input.height());
  
  for(int i = 0; i < input.height();i++){
	  for(int j = 0; j < input.width();j++){
		image[i*input.width() + j][0] = lab(j,i,0);  
		image[i*input.width() + j][1] = lab(j,i,1);  
		image[i*input.width() + j][2] = lab(j,i,2);  
	  }
  }
  int width= input.width();
  int height= input.height();
  
  int widthHolder = width;
  int heightHolder = height;
  
  
  
  for(int i = 0; i < widthHolder - atoi(argv[3]);i++){
	   seamer(image,height,width);
	   width--;
  }
  
  auto oldIm = image;
  //transpose the image
  for(int i = 0; i < height;i++){
	  for(int j = 0; j < width;j++){
		image[j*height+i] = oldIm[i*width+j];
	  }
  }
  swap(width,height);
  
  for(int i = 0; i < heightHolder - atoi(argv[4]);i++){
		seamer(image,height,width);
		width--;
  }
 
	oldIm = image;
  //transpose the image
  for(int i = 0; i < height;i++){
	  for(int j = 0; j < width;j++){
		image[j*height+i] = oldIm[i*width+j];
	  }
  }
  swap(width,height);
  
  CImg<double> output(height, width, input.depth(), input.spectrum(), 0);
  
  //CImg<double> output(atoi(argv[3]), atoi(argv[4]), input.depth(), input.spectrum(), 0);
  for (unsigned int i=0; i<output.height(); i++) {
	for (unsigned int j=0; j<output.width(); j++) {
	  output(j, i, 0) = image[i*output.width()+j][0];
	  output(j, i, 1) = image[i*output.width()+j][1];
	  output(j, i, 2) = image[i*output.width()+j][2];
	  
	  
	}
  }
  
  
  
  /*
  //display the energy
  CImg<double> output(width, height, input.depth(), input.spectrum(), 0);
  for (unsigned int i=0; i<output.height(); i++) {
	for (unsigned int j=0; j<output.width(); j++) {
	  double holder = pow(energies[i*output.width()+j]/MaxEnergy, 1.0/3.0) * 255.0;
	  output(i, j, 0) = holder;
	  output(i, j, 1) = holder;
	  output(i, j, 2) = holder;
	  
	}
  }
  */
  /*
  //display the cost
  double MaxCost = 0.0;
  
  for(int i = 0; i < costs.size(); i++){
	if(costs[i] > MaxCost){
		MaxCost = costs[i];
	}
  
  }
  
  
  CImg<double> output(width, height, input.depth(), input.spectrum(), 0);
  for (unsigned int i=0; i<output.height(); i++) {
	for (unsigned int j=0; j<output.width(); j++) {
	  double holder = pow(costs[i*output.width()+j]/MaxCost, 1.0/3.0) * 255.0;
	  output(i, j, 0) = holder;
	  output(i, j, 1) = holder;
	  output(i, j, 2) = holder;
	  
	}
  }
  */
  
  
 /*
  //display the seam
  CImg<double> output(width, height, input.depth(), input.spectrum(), 0);
  for (unsigned int i=0; i<output.height(); i++) {
	for (unsigned int j=0; j<output.width(); j++) {
	  output(j, i, 0) = 0.0;
	  output(j, i, 1) = 0.0;
	  output(j, i, 2) = 0.0;
		
	  if(j == seam[i]){
	  Eigen::Vector3d color = {150,100,200};
	  output(j, i, 0) = color[0];
	  output(j, i, 1) = color[1];
	  output(j, i, 2) = color[2];
	  }
	}
  }
  
  */
  
  /*
  CImg<double> output(atoi(argv[3]), atoi(argv[4]), input.depth(), input.spectrum(), 0);
  for (unsigned int i=0; i<output.width(); i++) {
	for (unsigned int j=0; j<output.height(); j++) {
	  output(i, j, 0) = image[i*output.height()+j][0];
	  output(i, j, 1) = image[i*output.height()+j][1];
	  output(i, j, 2) = image[i*output.height()+j][2];
	}
  }*/
  
  
  CImg<double> rgb = output.LabtoRGB();
 // CImg<double> rgb = output; //for debugging
  if (strstr(argv[2], "png"))
	rgb.save_png(argv[2]);
  else if (strstr(argv[2], "jpg"))
	rgb.save_jpeg(argv[2]);
  else if (strstr(argv[2], "bmp"))
	rgb.save_jpeg(argv[2]);

  //delete [] image;
  return 0;
}




