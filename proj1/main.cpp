#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <cmath>


using namespace std;


//write to ppm file
struct ppm{

  unsigned char*** pixels;

  void ppmSave(long &height, long &width){
    std::fstream fs;
    
    fs.open("hide.ppm", std::ios::out);

    //this wrotes the initial things for the ppm and then goes through each pixel with the right color
    fs << "P6" << std::endl;
    fs << width << " ";
    fs << height << std::endl;
    fs << 255 <<std::endl;;
    
     for(int i= 0; i < height; i++){

      for(int j = 0; j < width; j++){
	//1 = 1 byte from pixels location
	fs.write((char*)&pixels[i][j][0], 1);
	fs.write((char*)&pixels[i][j][1], 1);
	fs.write((char*)&pixels[i][j][2], 1);
      }
     }
    fs.close();

  }
  ppm(long &height, long &width){

    //collums
    pixels = new unsigned char**[height];

    for(int i= 0; i < height; i++){

      //rows
      pixels[i]= new unsigned char*[width];
      
      for(int j = 0; j < width; j++){

	//color (rgb)
	pixels[i][j] = new unsigned char[3];

	 //color (rgb) I made the default blue
          pixels[i][j][0] = 0;
          pixels[i][j][1] = 0;
          pixels[i][j][2] = 255;


	

      }

    }

  }

};


//creates the ray
class Ray{
public:
  Eigen::Vector3d origin, direction;
  Ray(const Eigen::Vector3d &o, const Eigen::Vector3d &d) : origin(o), direction(d){};


  
};

//used to keep track of hits, t, etc
class HitRecord{
public:
  bool hit = false;
  double t;
  Eigen::Vector3d color;

};



//intersection class to calculate the intersect

HitRecord intersect(const Ray &r, double t0, double t1, const std::vector<Eigen::Vector3d> &polys, Eigen::Vector3d &fil){
  Eigen::Matrix3d m;
  Eigen::Matrix3d a;
 
 Eigen::Matrix3d acopy;
	
 
  HitRecord hr;
  //t0 = 0 and t1 = closest known t value
  double alpha, betta, gama, t;

  double sol;
  //solution
   Eigen::Vector3d vex, solution;

   Eigen::Vector3d point;
   
/* this is how I came up with the math
	at first, i tried only following what i saw in class
	that didnt work so i followed the book
	it doesnt look pretty but i thought i should leave in the commented code 
	so you can see my process
*/
  
   //  m << polys[1] - polys[0], polys[2] - polys[0], -r.direction;

	 a << polys[0] - polys[1], polys[0] - polys[2], r.direction;
	 
	 sol = a.determinant();
	 
	 acopy = a;
     
	 solution = r.origin - polys[0];
     vex = m.inverse() * solution;

	acopy << polys[0] - r.origin, polys[0] - polys[2], r.direction;;
	
	
	betta = acopy.determinant()/sol;
	
	acopy = a;
	
	acopy <<  polys[0] - polys[1], polys[0] - r.origin, r.direction;
	
	gama = acopy.determinant()/sol;
	
	acopy = a;
	
	acopy<< polys[0] - polys[1], polys[0] - polys[2], polys[0] - r.origin;
	
	t = acopy.determinant()/sol;
	
	//cout << "b = " << hr.color[0] << " " << hr.color[1] << " "<< hr.color[2] << endl;
	
     //betta = vex[0];
     //gama = vex[1];
   //intesection value (compared with t0 and t1
     //t = vex[2];

     //alpha = 1 - betta - gama;


     //cout<< "beta " << betta << endl;
     // cout<< "gamma " << gama << endl;

     //     cout<< "m " << m << endl;
     // cout<< "polys " << polys[0] << " "<< polys[1] << " " << polys[2]<< endl;

  //----------------------------------intersection calculations-------------------------------------------------------------------------


       //check if NOT inside a triangle

     // if(0 < betta && 0 < gama && gama < 1 && (gama + betta) < 1  && t > t0 && t < t1){
/*
	if(betta > 0 && gama > 0 && t > 0 && betta + gama < 1){
		//inside the triangle
       point = polys[0] + betta * (polys[1] - polys[0]) + gama * (polys[2] - polys[0]);

       hr.color = fil;
       hr.t = t;
	   hr.hit = true;

       return true;
		
	}
	*/	 
	
	//checks from the book
	//this will keep track of if it hit and the t values
	//if it fullfills any of these if statements, then it did not intersect
     if(t < t0 || t > t1){
		
       hr.t = t;
	   hr.hit = false;
		return hr; 

				
	 }
	  if(gama < 0.0 || gama > 1.0){
		 
       hr.t = t;
	   hr.hit = false;
		return hr;   
	 }
	 
	if(betta < 0.0 || betta > (1.0 - gama)){
		
       hr.t = t;
	   hr.hit = false;
		return hr;   
		   
	   }
	     //pixel[betta][gama] = poly[0] + betta * (poly[1] + poly[0]) + gama * (poly[2] - poly[0]);

	     /*	       point = polys[0] + betta * (polys[1] - polys[0]) + gama * (polys[2] - polys[0]);

	       hr.color = fil;
	       hr.t = t;
	     */
	       //	       cout<<"point " << point<< endl;
	       //	       cout<<"hr.color " << hr.color<< endl;
	       //cout<< "hr.t " << hr.t<< endl;
	       //the ray hit
	       //if(r == point){
	       
	       //change color to red
	       //HitRecord();

	       //}
	       //}
	       //	       return true;
	//       return false;
	       //	   }
/*
	 }
       }
     }
  */     //inside the triangle
       //point = polys[0] + betta * (polys[1] - polys[0]) + gama * (polys[2] - polys[0]);

      
       hr.t = t;
	   hr.hit = true;

      
       return hr;
       
};





int main(int argc, char *argv[])
{

  cout<<"Working......"<< endl;
  
  
  
  //setting up everything to be used later on
  Eigen::Vector3d rgb, from, att, up, res, fil, vector, sphere, light;

  
  std::vector<std::vector<Eigen::Vector3d> > polys;

  
  double m,n,o, l, r, t, b;

  Eigen::Vector3d u,v,w;

  double angle, hither;


  long height, width;


  
  string line;
  //declare input stream and store the file
  ifstream nffInfo(argv[1]);


  //check if the file stored
  if (!nffInfo) {
    printf("Unable to open file %s", argv[1]);
    return 0;
  }


  string str;


  string::size_type sz;

 
    // Parse the nff file
  while(getline(nffInfo, line)){

     stringstream ss(line);

      char ch;
      //Background color.  A color is simply RGB with values between 0 and 1:
      if(line[0] == 'b'){

	ss>>ch>>rgb[0]>>rgb[1]>>rgb[2];
	//	cout << "b = " << rgb[0] << " " << rgb[1] << " "<< rgb[2] << endl;
      }

      //-----------------Viewpoint location.-------------------------------------------
      if(line[0] == 'v'){
	//	cout<<"v"<< endl;

	getline(nffInfo,line);
	stringstream jj(line);
      //the eye location in XYZ
	if(line[0] == 'f'){
	  jj>>str>>from[0]>>from[1]>>from[2];
	  //	  cout << "from = " << from[0] << " " << from[1] << " "<< from[2] << endl;

	}
	getline(nffInfo,line);
	stringstream kk(line);
	//a position to be at the center of the image, in XYZ world coordinates.  A.k.a. "lookat".
	if(line[0] == 'a'){

	  kk>>str>>att[0]>>att[1]>>att[2];
	  //cout << "at = " << att[0] << " " << att[1] << " "<< att[2] << endl;
	}

	getline(nffInfo,line);
	stringstream ll(line);
	//a vector defining which direction is up, as an XYZ vector.
	if(line[0] == 'u'){

	  ll>>str>>up[0]>>up[1]>>up[2];
	  // cout << "up = " << up[0] << " " << up[1] << " "<< up[2] << endl;

	}

      }
      //in degrees, defined as from the center of top pixel row to bottom pixel row and left column to right column.
      if(line[0] == 'a'){
	ss>>str>>angle;
	//cout << "angle = " << angle << endl;
      }

      //distance of the hither [near] plane (if any) from the eye.  Mostly needed for hidden surface algorithms.
      if(line[0] == 'h'){
	ss>>str>>hither;
	//	cout << "hither = " << hither << endl;

      }
      //in pixels, in x and in y.
      if(line[0] == 'r'){
	ss>>str>>res[0]>>res[1];

	height = res[1];
	width = res[0];
	//cout << "resolution = " << res[0] << " " << res[1] << endl;
      }

      //Positional light location
      if(line[0] == 'l'){
	ss>>str>>light[0]>>light[1]>>light[2];
	//cout << "l = " << light[0] << " " << light[1] << " " << light[2] << endl;
      }

      //fill color
      if(line[0] == 'f'){

	ss>>ch>>fil[0]>>fil[1]>>fil[2];
	//cout << "f = " << fil[0] << " " << fil[1] << " "<< fil[2] <<endl;//<< fil[0].Kd << " " << fil[0].Ks << " "<< fil[0].Shine << fil[0].T << " " << fil[0].iof << endl;

      }


      
      //polygon
      if(line[0] == 'p'){

	std::vector<Eigen::Vector3d> poly;
	//polycount++;
	int count;

	ss>>ch>>count;

	//cout << "p " << count << endl;


	for(int i = 0; i < count;i++){

	   getline(nffInfo,line);
	   stringstream jj(line);

	   jj>>vector[0]>>vector[1]>>vector[2];
	   poly.push_back(vector);
	   
	   //
	   //cout << vector[0] << " " << vector[1] << " " << vector[2] << endl;

		//cout << poly[0] << "\n" << poly[1] << "\n" << poly[2] << endl;
	   //   cout<<"vector " << vector<< endl;
	}


	polys.push_back(poly);
	//cout << polys[0][0]<< polys[0][1] << polys[0][2] <<endl;
	//cout<<"POLYS " << polys<< endl;
      }

      //sphere
      if(line[0] == 's'){

	ss>>ch>>sphere[0]>>sphere[1]>>sphere[2];//>>sphere[3];
	//cout << "s = " << sphere[0] << " " << sphere[1] << " "<< sphere[2] << endl;//" " << sphere[3] << endl;

      }

      //SKIP LINES
      if(line[0] == 'c'){

	continue;
      }

      //fill color
      if(line[0] == 'p' && line[1] == 'p'){
	continue;

      }


  }

//----------------------------- calulate ray from eye point through pixel and into the scene------------------------------------------
  unsigned char pixels[height][width][3];

//camera coordinates
  w = (from - att).normalized();
  u = up.cross(w).normalized();
  v = w.cross(u);//.normalized();

  

  //half of the images width
  double im_half_width = tan(M_PI*((double)angle/180.0)/2.0);
  double w_pixel_width = (im_half_width * 2.0)/(double)res[0];// * 2.0;
  double pixel_width = 2*tan(1.5*angle*M_PI/180.0)/width;
  double aspect = (double)res[0]/(double)res[1];

  
  
  //w_ is for all calculations for the world
  Eigen::Vector3d direction;
  Eigen::Vector3d w_direction;
  Eigen::Vector3d pixel;

  ppm output(height, width);

	int redcount = 0;
	int bluecount = 0; 
//coordinates on the image
	double imagex, imagey;
	
	Eigen::Vector3d pixelTarget;
	
	imagex = -0.5 * (width + 1) * pixel_width;
	imagey = 0.5 * (height + 1) * pixel_width;
  
  //goes through the y collum and then move over to the left for the x
  for(int i = 0; i < res[0]; i++){
    
	imagey = 0.5 * (height + 1) * pixel_width;
	
    for(int j = 0; j < res[1]; j++){
		
      l = -im_half_width;
      r = im_half_width;
      t = im_half_width / aspect;
      b = -im_half_width/ aspect;
      //the u and v of the pixel in the world space
      m = l + (w_pixel_width / 2.0) + j * w_pixel_width;
      n = t + (w_pixel_width / 2.0) - i * w_pixel_width;
      //d which is the 3rd coordinate in the camera system
      o = 1;
      
	  
	  //calculating pixel target
	  pixelTarget = imagex * u + imagey * v - w;
	  
	  direction = (from - pixelTarget).normalized();
      w_direction = m * u + n * v + o * -w;
      
      //pixel = from + direction;

      
      Ray r(from, direction);
	
	  Ray w_r(from, w_direction);
	  
      
      double min_t = numeric_limits<double>::infinity();
	  
	    HitRecord hr;
		
	
		hr.color = rgb;
		hr.hit = false;
		
	   //loop through all surface
      for(int k = 0; k < polys.size(); k++){
		  
		  
			auto result = intersect(w_r, 0.0, min_t, polys[k], fil);
			hr.t = result.t;
			hr.hit = result.hit;
		
			if(hr.t < min_t && hr.hit == true){
			   min_t = hr.t;
			   //pixel = hr.color;
			   hr.color = fil;
			 
			}
		 
	  

	   
		}	
		
		//assign the color and send to the ppm
		
	  output.pixels[i][j][0] = (unsigned char)(255 * hr.color(0));
	  output.pixels[i][j][1] = (unsigned char)(255 * hr.color(1));
	  output.pixels[i][j][2] = (unsigned char)(255 * hr.color(2));
	
		imagey -= pixel_width;

      }

		imagex += pixel_width;
	   
	  
    }

//creates the ppm and tells it how bug to be
    output.ppmSave(height, width);

	 
  return 0;
}
