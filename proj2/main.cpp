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

inline double sqr(double x) {return x*x;}
inline Eigen::Vector3d cross(const Eigen::Vector3d &x, const Eigen::Vector3d &y) {return x.cross(y);}
inline double dot(const Eigen::Vector3d &x, const Eigen::Vector3d &y) {return x.dot(y);}
inline double cross2d(const Eigen::Vector2d &x, const Eigen::Vector2d &y) {return x[0]*y[1]-x[1]*y[0];}

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

//handle fill from file
class Fill {
public: 
  Eigen::Vector3d color;
  double kd, ks, shine, transmittance, ior;
};


//handle light from the file
class Light {
public:
  Eigen::Vector3d p, c;
};



//creates the ray
class Ray{
public:
  Eigen::Vector3d origin, direction;
  double depth;
  Ray(const Eigen::Vector3d &o, const Eigen::Vector3d &d, double &de) : origin(o), direction(d), depth(de){};


  
};


//used to keep track of hits, t, etc
class HitRecord{
public:
  bool hit = false;
  double t, alpha, beta, gama;
  Fill f;
  Eigen::Vector3d color, p, n, v;

};

Eigen::Vector3d shade(const HitRecord &hr){
  return hr.f.color;
}

//intersection class to calculate the intersect

HitRecord intersect(const Ray &r, double t0, double t1, const std::vector<Eigen::Vector3d> &polys, Fill &fil){
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
	 i followed the book
	
*/
 // cout<<polys.size()<<endl;
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
	
	

  //----------------------------------intersection calculations-------------------------------------------------------------------------


       //check if NOT inside a triangle
	
	//checks from the book
	//this will keep track of if it hit and the t values
	//if it fullfills any of these if statements, then it did not intersect
     if(t < t0 || t > t1){
		
       hr.t = t;
	   hr.hit = false;
	   hr.f = fil;
		return hr; 

				
	 }
	  if(gama < 0.0 || gama > 1.0){
		 
       hr.t = t;
	   hr.hit = false;
	   hr.f = fil;
		return hr;   
	 }
	 
	if(betta < 0.0 || betta > (1.0 - gama)){
		
       hr.t = t;
	   hr.hit = false;
	   hr.f = fil;
		return hr;   
		   
	   }
	   
     //inside the triangle
       
	   point = polys[0] + betta * (polys[1] - polys[0]) + gama * (polys[2] - polys[0]);

      
       hr.t = t;
	   hr.hit = true;
	   hr.f = fil;
      
       return hr;
       
};



HitRecord SphereIntersect(const Ray &r, double t0, double t1, HitRecord &hr, Eigen::Vector3d &sphere, double &rad) {
  double ddotemc = dot(r.direction, r.origin-sphere);
  double d2 = r.direction.squaredNorm();

  double disc = sqr(ddotemc) - d2 * ((r.origin-sphere).squaredNorm() - rad*rad);

  if (disc < 0){
	  hr.hit = false;
	  return hr;
  }
  double root1 = (-ddotemc + sqrt(disc)) / d2;
  double root2 = (-ddotemc - sqrt(disc)) / d2;

  double t = root1;
  if (root1 < 0 || (root2 > 0 && root2 < root1)) t = root2;
  if (t < t0 || t > t1){
	  hr.hit = false;
	  return hr;
  }
  
  hr.t = t;
  hr.p = r.origin + t * r.direction;
  hr.n = (hr.p - sphere) / rad;
  hr.hit = true;
  return hr;
}





//--------taken from progs .h

class Surface {
public:
  //virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const = 0;
  //virtual Eigen::Vector3d normal(const Eigen::Vector3d &x) const = 0;
  virtual ~Surface() {};
};

class Triangle : public Surface {
protected:
  Eigen::Vector3d a,b,c;
public:
  Triangle(const Eigen::Vector3d &_a, const Eigen::Vector3d &_b, const Eigen::Vector3d &_c) : a(_a), b(_b), c(_c) {};
  //virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
  // virtual Eigen::Vector3d normal(const Eigen::Vector3d &x) const;
};

class TrianglePatch : public Triangle {
  Eigen::Vector3d n1, n2, n3;
public:
  TrianglePatch(const Eigen::Vector3d &_a, const Eigen::Vector3d &_b, const Eigen::Vector3d &_c,
	  const Eigen::Vector3d &_n1, const Eigen::Vector3d &_n2, const Eigen::Vector3d &_n3) 
	: Triangle(_a,_b,_c), n1(_n1), n2(_n2), n3(_n3) {};
  // virtual Eigen::Vector3d normal(const Eigen::Vector3d &x) const;
};

class Poly : public Surface {
  std::vector<Eigen::Vector3d> verts;
public:
  Poly(const std::vector<Eigen::Vector3d> &_verts) : verts(_verts) {}; 
  //virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
  // virtual Eigen::Vector3d normal(const Eigen::Vector3d &x) const;
};

class PolyPatch : public Poly {
  std::vector<Eigen::Vector3d> normals;
public:
  PolyPatch(const std::vector<Eigen::Vector3d> &_verts, const std::vector<Eigen::Vector3d> &_normals) : Poly(_verts), normals(_normals) {}; 
};

class Sphere : public Surface {
  Eigen::Vector3d c;
  double rad;
public:
  Sphere(const Eigen::Vector3d &_c, double _r) : c(_c), rad(_r) {};
 // virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
  // virtual Eigen::Vector3d normal(const Eigen::Vector3d &x) const;
};


///-------------------

int main(int argc, char *argv[])
{

  cout<<"Working......"<< endl;
  
  
  
  std::vector<Fill> fills;
  
  
  //setting up everything to be used later on
  Eigen::Vector3d rgb, from, att, up, res, vector, nvec, sphere, light;
  
  std::vector<std::vector<Eigen::Vector3d> > polys; 
  std::vector<std::pair<Surface *, Fill> > polycount;

  std::vector<std::pair<Surface *, Fill> > surfaces;
  std::vector<Light> lights;
  
  double m,n,o, l, r, t, b;

  Eigen::Vector3d u,v,w;

  double angle, hither;

  long height, width;
  
  bool coloredlights = false;
  double aperture = 0.0;
  int samples = 1;
  int maxraydepth = 5;
  bool color = false;
  bool phong = false;
  bool stratified = false;
  bool reflections = true;
  bool shadows = true;
  bool verbose = false;
  int mode = 3;
  double shadowbias; 
  
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

	Fill fil;
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
		  
		Light l;
		ss>>ch>>l.p[0]>>l.p[1]>>l.p[2];
		
		//check if there are colored lights
		if (!ss.eof()) {
		ss>>l.c[0]>>l.c[1]>>l.c[2];
		coloredlights = true;
	  }
		
		lights.push_back(l);
	//ss>>str>>light[0]>>light[1]>>light[2];
	//cout << "l = " << light[0] << " " << light[1] << " " << light[2] << endl;
      }
	//if there are NOT colored lights  
	  if (!coloredlights) {
		for (unsigned int i=0; i<lights.size(); i++) {
		  lights[i].c[0] = 1.0/sqrt(lights.size());
		  lights[i].c[1] = 1.0/sqrt(lights.size());
		  lights[i].c[2] = 1.0/sqrt(lights.size());
		}
	  }
	  
	  shadowbias = 1e-6;
	  samples = 1;
	  aperture = 0.0;

      //fill color
      if(line[0] == 'f'){
       ss>>ch>>fil.color[0]>>fil.color[1]>>fil.color[2]>>fil.kd>>fil.ks>>fil.shine>>fil.transmittance>>fil.ior;
		fills.push_back(fil);
	//ss>>ch>>fil[0]>>fil[1]>>fil[2];
	//cout << "f = " << fil.color[0] << " " << fil.color[1] << " "<< fil.color[2] << " " << fil.kd << " " << fil.ks << " "<< fil.shine << fil.transmittance << " " << fil.ior << endl;

      }

		bool patch = false;
      
      //polygon
      if(line[0] == 'p'){
		  
		  if(line[1] == 'p'){
			patch = true;
			ss>>ch;
		  }
		std::vector<Eigen::Vector3d> poly; //vertices
		std::vector<Eigen::Vector3d> normals;
		
		
		int count; //nverts

		ss>>ch>>count;


	for(int i = 0; i < count;i++){

	   getline(nffInfo,line);
	   stringstream jj(line);

		
		if(patch)jj>>vector[0]>>vector[1]>>vector[2]>>nvec[0]>>nvec[1]>>nvec[2];
		   
		
		else jj>>vector[0]>>vector[1]>>vector[2];
		
		poly.push_back(vector);
		
		if (patch) {
			  nvec.normalize();
			  normals.push_back(nvec);
			}
	
		}

//-------------------------

		bool makeTriangles = false;
	  if (poly.size() == 3) {
	        makeTriangles = true;
		if (patch) {
		  surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(poly[0], poly[1], poly[2], 
					  normals[0], normals[1], normals[2]), fills.back()));
		  polys.push_back(poly);
		  
		} else {
		  surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(poly[0], poly[1], poly[2]), fills.back()));
		 polys.push_back(poly);
		
		}
		//poly has 4 sides
	  } else if (poly.size() == 4) {
		 		
		Eigen::Vector3d n0 = cross(poly[1] - poly[0], poly[2] - poly[0]);
		Eigen::Vector3d n1 = cross(poly[2] - poly[1], poly[3] - poly[1]);
		Eigen::Vector3d n2 = cross(poly[3] - poly[2], poly[0] - poly[2]);
		Eigen::Vector3d n3 = cross(poly[0] - poly[3], poly[1] - poly[3]);
		
		if (dot(n0,n1) > 0 && dot(n0,n2) > 0 && dot(n0,n3) > 0) {
			
		  makeTriangles = true;
		  if (patch) {
			 // cout<<"in if" << endl;
			surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(poly[0], poly[1], poly[2], 
						normals[0], normals[1], normals[2]), fills.back()));
			surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(poly[0], poly[2], poly[3], 
						normals[0], normals[2], normals[3]), fills.back()));
						
			polys.push_back(poly);
			polys.push_back(poly);
			
		  } else {
			 // cout<<"else " << endl;
			surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(poly[0], poly[1], poly[2]), fills.back()));
			surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(poly[0], poly[2], poly[3]), fills.back()));			
			
			polys.push_back(poly);
			polys.push_back(poly);
			
		  }
		}
	  }
	  else{
		  //fan them
		  //auto *first_vertex = poly[0];
		  for(int i = 0; i < poly.size() - 1; i++){
			  surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(poly[0], poly[1], poly[i+1]), fills.back()));
			  polys.push_back(poly);
			  
		  }
		}	
		 if (!makeTriangles) {
			if (patch) {
			  surfaces.push_back(std::pair<Surface *, Fill>(new PolyPatch(poly, normals), fills.back()));
			  polys.push_back(poly);
			  
			} else {
			  surfaces.push_back(std::pair<Surface *, Fill>(new Poly(poly), fills.back()));
			  polys.push_back(poly);
			 
			} 
		  } 
	 
//------------------
		
		
		//polys.push_back(poly);
		//cout << polys[0][0]<< polys[0][1] << polys[0][2] <<endl;
		//cout<<"POLYS " << polys<< endl;
      }

      //sphere
      if(line[0] == 's'){
		std::vector<Eigen::Vector3d> spheres;

		double r;
		ss>>ch>>sphere[0]>>sphere[1]>>sphere[2]>>r;
		surfaces.push_back(std::pair<Surface *, Fill>(new Sphere(sphere,r), fil));
		spheres.push_back(sphere);
		polys.push_back(spheres);
	//cout << "s = " << sphere[0] << " " << sphere[1] << " "<< sphere[2] << endl;//" " << sphere[3] << endl;

      }

      //SKIP LINES
      if(line[0] == 'c'){

	continue;
      }
/*
      //fill color
      if(line[0] == 'p' && line[1] == 'p'){
	continue;

      }
*/

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
  
  cout<<polys.size()<<endl;
  cout<<surfaces.size()<<endl;
  
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

      
      Ray r(from, direction, o);
	
	  Ray w_r(from, w_direction, o);
	  
      
      double min_t = numeric_limits<double>::infinity();
	  
	    HitRecord hr;
		
	
		hr.color = rgb;
		hr.hit = false;
		
	   //loop through all surface
      for(int k = 0; k < polys.size(); k++){
		  
		  
			auto result = intersect(w_r, 0.0, min_t, polys[k], surfaces[k].second);
			hr.t = result.t;
			hr.hit = result.hit;
		
			if(hr.t < min_t && hr.hit == true){
			   min_t = hr.t;
			   //pixel = hr.color;
			   hr.color = surfaces[k].second.color;
			   
			 
			}
		 	  shade(hr);
			if (hr.hit) {
				color = shade;
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
