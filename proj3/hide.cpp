#include "trace.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <getopt.h>
#define cimg_display 0
#include "CImg.h"
#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else
#include <values.h>
#define MAX DBL_MAX
#endif
 
using namespace cimg_library;
using namespace std;

Fragment::Fragment(double a, Eigen::Vector3d b){
	z = a;
	color = b;
	
} 


// return the determinant of the matrix with columns a, b, c.
double det(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c) {
  return a[0]* (b[1] * c[2] - c[1] * b[2]) +
	b[0] * (c[1] * a[2] - a[1] * c[2]) +
	c[0] * (a[1] * b[2] - b[1] * a[2]);
}



inline double sqr(double x) {return x*x;}
inline Eigen::Vector3d cross(const Eigen::Vector3d &x, const Eigen::Vector3d &y) {return x.cross(y);}
inline double dot(const Eigen::Vector3d &x, const Eigen::Vector3d &y) {return x.dot(y);}
inline double cross2d(const Eigen::Vector2d &x, const Eigen::Vector2d &y) {return x[0]*y[1]-x[1]*y[0];}

inline double area(const Eigen::Vector2d &a, const Eigen::Vector2d &b, const Eigen::Vector2d &c) {
	Eigen::Vector2d p1 = b - c;
	double area = cross2d(p1, c-a);
	return area/2;}


bool Triangle::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
  Eigen::Vector3d ba = a-b;
  Eigen::Vector3d ca = a-c;
  Eigen::Vector3d ea = a-r.e;
  double detA = det(ba, ca, r.d);
  double beta = det(ea, ca, r.d)/detA;
  if (beta < 0 || beta > 1) return false;
  double gamma = det(ba, ea, r.d)/detA;
  if (gamma < 0.0 || gamma > 1.0-beta) return false;
  double t = det(ba, ca, ea)/detA;
  if (t < t0 || t > t1) return false;
  hr.t = t;
  hr.p = r.e + t * r.d;
  hr.n = cross(ba,ca);
  hr.n.normalize();
  hr.alpha = 1.0 - beta - gamma;
  hr.beta = beta;
  hr.gamma = gamma;
  return true;
}


bool Sphere::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
  double ddotemc = dot(r.d, r.e-c);
  double d2 = r.d.squaredNorm();

  double disc = sqr(ddotemc) - d2 * ((r.e-c).squaredNorm() - rad*rad);

  if (disc < 0) return false;
  double root1 = (-ddotemc + sqrt(disc)) / d2;
  double root2 = (-ddotemc - sqrt(disc)) / d2;

  double t = root1;
  if (root1 < 0 || (root2 > 0 && root2 < root1)) t = root2;
  if (t < t0 || t > t1) return false;
  
  hr.t = t;
  hr.p = r.e + t * r.d;
  hr.n = (hr.p - c) / rad;
  return true;
}

Eigen::Vector2d project(const Eigen::Vector3d &x, int projectDir) {
  switch (projectDir) {
  case 0:
    return Eigen::Vector2d(x[1],x[2]);
  case 1:
    return Eigen::Vector2d(x[0],x[2]);
  case 2:
    return Eigen::Vector2d(x[0],x[1]);
  }
  return Eigen::Vector2d(1.0, 1.0);
}

bool Poly::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
  int projectDir;
  Eigen::Vector3d n = cross(verts[1]-verts[0], verts[2]-verts[0]);
  n.normalize();

  double t = -(dot(r.e, n) - dot(verts[0], n)) / dot(r.d,n);
  if (t < t0 || t > t1) return false;

  Eigen::Vector3d p = r.e + t*r.d;

  if (fabs(n[0]) > fabs(n[1]) && fabs(n[0]) > fabs(n[2])) {
	projectDir = 0;
  } else if (fabs(n[1]) > fabs(n[2])) {
	projectDir = 1;
  } else {
	projectDir = 2;
  }
  
  Eigen::Vector2d p2 = project(p, projectDir);

  Eigen::Vector2d bbMin = project(verts[0], projectDir);
  Eigen::Vector2d bbMax = bbMin;
  for (unsigned int i=1; i<verts.size(); i++) {
	Eigen::Vector2d v = project(verts[i], projectDir);
	if (v[0] < bbMin[0]) bbMin[0] = v[0];
	if (v[0] > bbMax[0]) bbMax[0] = v[0];
	if (v[1] < bbMin[1]) bbMin[1] = v[1];
	if (v[1] > bbMax[1]) bbMax[1] = v[1];
  }
  
  if (p2[0] < bbMin[0]) return false;
  if (p2[1] < bbMin[1]) return false;
  if (p2[0] > bbMax[0]) return false;
  if (p2[1] > bbMax[1]) return false;

  Eigen::Vector2d dir(sqrt(2), sqrt(2));
  int count = 0;
  for (unsigned int i=0; i<verts.size(); i++) {
	Eigen::Vector2d a = project(verts[i], projectDir);
	Eigen::Vector2d b = project(verts[(i+1) % verts.size()], projectDir);
	Eigen::Vector2d ab = b-a;
	double t2 = cross2d(a-p2, ab / cross2d(dir, ab));
	if (t2 < 0.0) continue;
	double alpha = cross2d(a-p2, dir / cross2d(dir, ab));
	if (alpha > 0.0 && alpha < 1.0) count++;
  }

  if (count % 2 == 0) {
	return false;
  }

  hr.t = t;
  hr.p = p;
  hr.n = n;
  return true;
}

Tracer::Tracer(const std::string &fname) {
  std::ifstream in(fname.c_str(), std::ios_base::in);
  std::string line;
  char ch;
  Fill fill;
  bool coloredlights = false;

  while (in) {
	getline(in, line);
	switch (line[0]) {
	case 'b': {
	  std::stringstream ss(line);
	  ss>>ch>>bcolor[0]>>bcolor[1]>>bcolor[2];
	  break;
	}

	case 'v': {
	  getline(in, line);
	  std::string junk;
	  std::stringstream fromss(line);
	  fromss>>junk>>eye[0]>>eye[1]>>eye[2];

	  getline(in, line);
	  std::stringstream atss(line);
	  atss>>junk>>at[0]>>at[1]>>at[2];

	  getline(in, line);
	  std::stringstream upss(line);
	  upss>>junk>>up[0]>>up[1]>>up[2];

	  getline(in, line);
	  std::stringstream angless(line);
	  angless>>junk>>angle;

	  getline(in, line);
	  std::stringstream hitherss(line);
	  hitherss>>junk>>hither;

	  getline(in, line);
	  std::stringstream resolutionss(line);
	  resolutionss>>junk>>res[0]>>res[1];
	  break;
	}

	case 'p': {
	  bool patch = false;
	  std::stringstream ssn(line);
	  unsigned int nverts;
	  if (line[1] == 'p') {
		patch = true;
		ssn>>ch;
	  }
	  ssn>>ch>>nverts;
	  std::vector<Eigen::Vector3d> vertices;
	  std::vector<Eigen::Vector3d> normals;
	  for (unsigned int i=0; i<nverts; i++) {
		getline(in, line);
		std::stringstream ss(line);
		Eigen::Vector3d v,n;
		if (patch) ss>>v[0]>>v[1]>>v[2]>>n[0]>>n[1]>>n[2];
		else ss>>v[0]>>v[1]>>v[2];
		vertices.push_back(v);
		if (patch) {
		  n.normalize();
		  normals.push_back(n);
		}
	  }

	  bool makeTriangles = false;
	  if (vertices.size() == 3) {
	        makeTriangles = true;
		if (patch) {
		  surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2], 
					  normals[0], normals[1], normals[2]), fill));
		} else {
		  surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
		}
	  } else if (vertices.size() == 4) {
		Eigen::Vector3d n0 = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
		Eigen::Vector3d n1 = cross(vertices[2] - vertices[1], vertices[3] - vertices[1]);
		Eigen::Vector3d n2 = cross(vertices[3] - vertices[2], vertices[0] - vertices[2]);
		Eigen::Vector3d n3 = cross(vertices[0] - vertices[3], vertices[1] - vertices[3]);
		if (dot(n0,n1) > 0 && dot(n0,n2) > 0 && dot(n0,n3) > 0) {
		  makeTriangles = true;
		  if (patch) {
			surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2], 
						normals[0], normals[1], normals[2]), fill));
			surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[2], vertices[3], 
						normals[0], normals[2], normals[3]), fill));
		  } else {
			surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
			surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[2], vertices[3]), fill));
		  }
		}
	  }
	  if (!makeTriangles) {
	    if (patch) {
	      surfaces.push_back(std::pair<Surface *, Fill>(new PolyPatch(vertices, normals), fill));
	    } else {
	      surfaces.push_back(std::pair<Surface *, Fill>(new Poly(vertices), fill));
	    }
	  }
	  break;
	}

	case 's' : {
	  std::stringstream ss(line);
	  Eigen::Vector3d c;
	  double r;
	  ss>>ch>>c[0]>>c[1]>>c[2]>>r;
	  surfaces.push_back(std::pair<Surface *, Fill>(new Sphere(c,r), fill));
	  break;
	}
	  
	case 'f' : {
	  std::stringstream ss(line);
	  ss>>ch>>fill.color[0]>>fill.color[1]>>fill.color[2]>>fill.kd>>fill.ks>>fill.shine>>fill.transmittance>>fill.ior;
	  break;
	}

	case 'l' : {
	  std::stringstream ss(line);
	  Light l;
	  ss>>ch>>l.p[0]>>l.p[1]>>l.p[2];
	  if (!ss.eof()) {
		ss>>l.c[0]>>l.c[1]>>l.c[2];
		coloredlights = true;
	  }
	  lights.push_back(l);
	  break;
	}

	default:
	  break;
	}
  }


  if (!coloredlights) {
    for (unsigned int i=0; i<lights.size(); i++) {
      lights[i].c[0] = 1.0/sqrt(lights.size());
      lights[i].c[1] = 1.0/sqrt(lights.size());
      lights[i].c[2] = 1.0/sqrt(lights.size());
    }
  }
  im = new Eigen::Vector3d[res[0]*res[1]];
  shadowbias = 1e-6;
  samples = 1;
  aperture = 0.0;
  

}

Tracer::~Tracer() {
  if (im) delete [] im;
  for (unsigned int i=0; i<surfaces.size(); i++) delete surfaces[i].first;
}

//------------------------------my code-------------------------------------------------
Eigen::Vector3d Tracer::shade(const HitRecord &hr,const Light &l, Fill &f,const Ray &ray, int bounceLimits) const {
	//all of the math is taken from the book or lecture
	HitRecord hr2;
	bool inShadow = false;
	bool hit = false;
	Eigen::Vector3d color = {0,0,0};
	Triangle *tri;
	//direction from the eye to the point
	Eigen::Vector3d d = (ray.e-hr.p).normalized();
	Eigen::Vector3d v = -d;
	//reflected vector
	Eigen::Vector3d r = d - 2 * (d.dot(hr.n)) * hr.n;
	r = r.normalized();
	
	double bias = 0.001;
	
	//is it in shadow?
	for(int i = 0; i < lights.size(); i++){
		//light direction
		Eigen::Vector3d ld = (lights[i].p-hr.p).normalized();
		inShadow = false;
		double t1 = (lights[i].p-hr.p).norm();
		//go through each surface and see if its in shadow
		for(int j = 0; j < surfaces.size(); j++){
			Ray shadowRay(hr.p, ld);
			const std::pair<Surface *, Fill> &s  = surfaces[j];
			//check to see if its in shadow
			if (s.first->intersect(shadowRay, bias, t1, hr2)) {
				inShadow = true;
			}
		}
		//color
		if(!inShadow){
			
		//diffuse and spectular
			double nl = hr.n.dot(ld);
			double vr = v.dot(r);
			
			if(nl < 0){
				nl = 0;  
			}
			if(vr < 0){
				vr = 0;
			}
			//from the book
			color[0] += (f.color[0] * f.kd * nl + f.ks*pow(vr, f.shine)) * lights[i].c[0];
			color[1] += (f.color[1] * f.kd * nl + f.ks*pow(vr, f.shine)) * lights[i].c[1];
			color[2] += (f.color[2] * f.kd * nl + f.ks*pow(vr, f.shine)) * lights[i].c[2];
		
			
		}
		
	}
	
	double t1 = MAX;
	
	//compute reflected colors
	if(f.ks > 0 && bounceLimits > 0){
		hit = false;
		Ray reflectRay(hr.p, -r);
		Eigen::Vector3d reflectColor(bcolor);
		//go through each surface to see if it needs to be reflected
		for(int j = 0; j < surfaces.size(); j++){
			const std::pair<Surface *, Fill> &s  = surfaces[j];
			//save all of the info
			if (s.first->intersect(reflectRay, bias, t1, hr2)) {
				t1 = hr2.t;
			    hr2.f = s.second;
			    hr2.raydepth = reflectRay.depth;
			    hr2.v = reflectRay.e - hr2.p;
			    hr2.v.normalize();
			    hit = true;
			}
		}
		//if it hit a surface we want it to reflect
		if(hit){
			//recursive call (math from the book)
			color += f.ks * this->shade(hr2, lights[0], hr2.f, reflectRay, bounceLimits - 1);			
		}
		else{
			//we didnt hit so nothing to reflect
			color += f.ks * reflectColor;
		}
		
		
	}
	
  return color;
 
}
//----------------------------------------------------------(i do change a few things in the header and in cstray but otherwise its the same)


//access a,b,c from triangle (vertices)
Eigen::Vector3d Triangle::getVertex(int i){
		switch(i){
				case 0:
				return a;
				break;
				case 1:
				return b;
				break;
				case 2:
				return c;
				break;
		}
	return (Eigen::Vector3d){0,0,0};
}

//access n1,n2,n3 from Triangle(normals)
Eigen::Vector3d Triangle::getNormal(int i){
		
	return ((b - a).cross(c - a)).normalized();
}

//access n1,n2,n3 from TrianglePatch(normals)
Eigen::Vector3d TrianglePatch::getNormal(int i){
		switch(i){
				case 0:
				return n1;
				break;
				case 1:
				return n2;
				break;
				case 2:
				return n3;
				break;
		}
	return (Eigen::Vector3d){0,0,0};
}


Eigen::Vector3d Tracer::castRay(const Ray &r, double t0, double t1) const {
  HitRecord hr;
  Eigen::Vector3d color(bcolor);
  
  bool hit = false;
  for (unsigned int k=0; k<surfaces.size(); k++) {
	const std::pair<Surface *, Fill> &s  = surfaces[k];
	if (s.first->intersect(r, t0, t1, hr)) {
	  t1 = hr.t;
	  hr.f = s.second;
	  hr.raydepth = r.depth;
	  hr.v = r.e - hr.p;
	  hr.v.normalize();
	  hit = true;
	}
  }

  if (hit) {
	color = this->shade(hr, lights[0], hr.f, r, 5);
  }
  return color;
}


//compares fragments
bool compFrag(Fragment a, Fragment b){
	return a.z<b.z;
	
}




void Tracer::createImage() {
	
  // set up coordinate system
  Eigen::Vector3d w = eye - at;
  w /= w.norm();
  Eigen::Vector3d u = cross(up,w);
  u.normalize();
  Eigen::Vector3d v = cross(w,u);

  Eigen::Vector3d color = {0,0,0};
  //std::cout<<u<<" "<<v<<" "<<w<<std::endl;
  
  double d = (eye - at).norm();
  double h = tan((M_PI/180.0) * (angle/2.0)) * d;
  double increment = (2*h) / res[0];
  double l = -h + 0.5*increment;
  double t = h*(((double)res[1])/res[0]) - 0.5*increment;

  if (verbose) {
    std::cout<<"u: ["<<u[0]<<", "<<u[1]<<", "<<u[2]<<"]"<<std::endl;
    std::cout<<"v: ["<<v[0]<<", "<<v[1]<<", "<<v[2]<<"]"<<std::endl;
    std::cout<<"w: ["<<w[0]<<", "<<w[1]<<", "<<w[2]<<"]"<<std::endl;
    std::cout<<"d: "<<d<<std::endl;
    std::cout<<"deltax / h: "<<increment<<std::endl;
    std::cout<<"left limit: "<<l<<std::endl;
    std::cout<<"top limit: "<<t<<std::endl;
  }

  

  /*for (unsigned int j=0; j<res[1]; j++) {
    for (unsigned int i=0; i<res[0]; i++, pixel++) {
      double x = l + i*increment;
      double y = t - j*increment;
      Eigen::Vector3d dir = x*u + y*v - d*w;
      (*pixel) = Eigen::Vector3d(0.0,0.0,0.0);
      for (int k=0; k<samples; k++) {
		for (int l=0; l<samples; l++) {
		  Eigen::Vector3d origin = eye;
		  Eigen::Vector3d imagept = eye + dir;
		  Ray r(origin, imagept-origin);
		  r.d.normalize();
		//  (*pixel) += castRay(r, hither, MAX) / (samples*samples);
		}
      }
    }
  }
  */
//----------------------------------proj3---------------------------------------
Frustum frustum;
double alpha,beta,gamma;
double aspect = res[1]/res[0];

double H = frustum.h = tan((M_PI/180.0) * (angle/2.0)) * hither;
double L = frustum.l = -H;
double R = frustum.r = H;
double B = frustum.b = -H/aspect;
double T = frustum.t = H/aspect;
double n = frustum.n = -hither;
double f = frustum.f = -1000.0*hither; 
//the vectors coordinates



//lower ceiling(min) and upper ceiling(max)
Eigen::Vector2d lc,uc;

//matrix multiply
Eigen::Matrix4d Mvp;
Eigen::Matrix4d Mper;
Eigen::Matrix4d Mcam;

Mvp<<res[0]/2, 0.0, 0.0, (res[0]-1.0)/2.0,
	0.0, res[1]/2.0,0.0, (res[1]-1.0)/2.0,
	0.0, 0.0, 1.0, 0.0,
	0.0, 0.0, 0.0, 1.0;

Mper<<(2*n)/(R-L), 0.0, (L+R)/(L-R), 0.0,
	0.0, (2*n)/(T-B),(B+T)/(B-T), 0.0,
	0.0, 0.0, (f+n)/(n-f), (2*f*n)/(f-n),
	0.0, 0.0, 1.0, 0.0;
	
Mcam<<u[0], v[0], w[0], eye[0],
	u[1], v[1], w[1], eye[1],
	u[2], v[2], w[2], eye[2],
	0.0, 0.0, 0.0, 1.0;
Mcam = Mcam.inverse().eval();

Eigen::Matrix4d M = Mvp * Mper * Mcam; //p153
	

Eigen::Vector3d V[3];
Eigen::Vector3d N[3];
Eigen::Vector4d v_prime[3];

Eigen::Vector3d col[3]; //colors

Eigen::Vector3d abgColor;
std::vector<Fragment> *fragments = new std::vector<Fragment>[res[0]*res[1]];


	for(int k = 0; k < surfaces.size(); k ++){
		//go through each vertex of the triangle
		
		for(int i = 0; i < 3; i++){
			V[i] = ((Triangle *) (surfaces[k].first))->getVertex(i);
			Eigen::Vector4d taco;
			//taco is a placeholder name
			taco[0] = V[i][0];
			taco[1] = V[i][1];
			taco[2] = V[i][2];
			taco[3] = 1.0;
			
			col[i] = {0,0,0};
			
			//N[i] = ((TrianglePatch *) (surfaces[k].first))->getNormal(i);
			N[i] = ((Triangle *) (surfaces[k].first))->getNormal(i);
			
			
			//shading part
			for(int g = 0; g < lights.size(); g++){
				double lightIntensity = 1/sqrt(lights.size());
				//for sepcular
				//direction from the eye to the point
				Eigen::Vector3d dir = (eye-V[i]).normalized();
				Eigen::Vector3d ve = -dir;
				//reflected vector
				Eigen::Vector3d ree = dir - 2 * (dir.dot(N[i])) * N[i];
				ree = ree.normalized();
				
				//spectular
				double vree = ve.dot(ree);
				
				if(vree < 0){
					vree = 0;
				}
				//from the book
				col[i][0] += (surfaces[k].second.ks*pow(vree, surfaces[k].second.shine)) * lightIntensity;
				col[i][1] += (surfaces[k].second.ks*pow(vree, surfaces[k].second.shine)) * lightIntensity;
				col[i][2] += (surfaces[k].second.ks*pow(vree, surfaces[k].second.shine)) * lightIntensity;
								
				Eigen::Vector3d lightDirection = lights[g].p - V[i];
				lightDirection.normalize();
				double LD = max(0.0, N[i].dot(lightDirection)); 
				col[i] += LD * surfaces[k].second.kd * lightIntensity * surfaces[k].second.color;
			}
			
			
			v_prime[i] = M * taco;
			
			v_prime[i] = v_prime[i] / v_prime[i][3]; //homegenous devide
		}

		lc[0] = v_prime[0][0];
		uc[0] = v_prime[0][0];

		//loop over v_prime and take floor and ceil for min and max
		//our boundingbox
		for(int i = 1; i<3; i++){
			lc[0] = std::min(lc[0], v_prime[i][0]);
			lc[1] = std::min(lc[1], v_prime[i][1]);
			uc[0] = std::max(uc[0], v_prime[i][0]);
			uc[1] = std::max(uc[1], v_prime[i][1]);
		}

		lc[0] = std::max(std::floor(lc[0]), 0.0);
		lc[1] = std::max(std::floor(lc[1]), 0.0);
		uc[0] = std::min(std::ceil(uc[0]), (double) res[0]);
		uc[1] = std::min(std::ceil(uc[1]), (double) res[1]);


		//vertex processing 


		//persperctive devide
		//compute colors


		Eigen::Vector2d p;
		Eigen::Vector2d vp[3];

		//rasterizer 
		//for(int i=(int)lc[0]; i <= (int)uc[0]; i++){
			//for(int j= (int)lc[1]; j<= (int)uc[1]; j++){
		for(int i = 0; i < res[0]; i++){
			for(int j = 0; j<res[1];j++){
			//for(int j= (int)lc[1]; j<= (int)uc[1]; j++){
				//for(int i=(int)lc[0]; i <= (int)uc[0]; i++){
				//position of in the picture
				p[0] = i; //* (res[0] + lc[1]) + j;
				p[1] = j;
				
				vp[0][0] = v_prime[0][0];
				vp[1][0] = v_prime[1][0];
				vp[1][1] = v_prime[1][1];
				vp[0][1] = v_prime[0][1];
				vp[2][0] = v_prime[2][0];
				vp[2][1] = v_prime[2][1];
				
				
				
				//compute alpha,beta,gamma see page 165
				alpha = cross2d(vp[1] - p,vp[2] - p)/cross2d(vp[1] - vp[0],vp[2] - vp[0]);
				beta = cross2d(vp[2] - p,vp[0] - p)/cross2d(vp[1] - vp[0],vp[2] - vp[0]);
				gamma = cross2d(vp[0] - p,vp[1] - p)/cross2d(vp[1] - vp[0],vp[2] - vp[0]);
				
				abgColor = alpha * col[0] + beta * col[1] + gamma * col[2];//surfaces[k].second.color;
				
				//if(i == 396 && j == 319){
					//	cout<<"alpha = "<< alpha<< "beta = " << beta<< "gamma = " << gamma<< "p = "<< p[0] << p[1]<< endl;
						//cout<<fragments[i+j*res[1]].size()<<endl;
				//}
				
				
				//check if its in the triangle
				if(alpha > 0 && beta > 0 && gamma > 0){
					
					//generate fragment
					fragments[(res[0] * res[1] - res[1]) - j*res[1] + i].push_back(Fragment(alpha*v_prime[0][2]+beta*v_prime[1][2]+gamma*v_prime[2][2], abgColor));
					//fragments[i+j*res[1]].push_back(Fragment(alpha*v_prime[0][2]+beta*v_prime[1][2]+gamma*v_prime[2][2], abgColor));
					
					//if(i == 396 && j == 319){
						//cout<<"alpha = "<< alpha<< "beta = " << beta<< "gamma = " << gamma<< "p = "<< p[0] << p[1]<< endl;
						//cout<<fragments[i+j*res[1]].size()<<endl;
						//}
					
				}
			}
		}
	}
		//fragment processing
		Eigen::Vector3d *pixel = im;
		//blending
		for(unsigned int j = 0; j < res[1]; j++){
			for(unsigned int i = 0; i < res[0]; i++, pixel++){
			
				unsigned int index = i+j*res[1]; 
				//unsigned int index = (res[0] * res[1] - res[1]) - j*res[1] + i;
				double minz = -DBL_MAX;
				color = bcolor;//background;
				
				//std::sort(fragments[index].begin(), fragments[index].end(), compFrag);
				
				//if(i == 396 && j == 319){
				//cout<<"fragments "<< fragments[index].size()<<endl;
				//}
				
				//loop through fragments
				for(int m = 0; m < fragments[index].size(); m++){
					
						if(fragments[index][m].z > minz){
							color= fragments[index][m].color;
							minz = fragments[index][m].z;
						}
					
				}
				(*pixel) = color;
			}	
		}

		//frame buffer image

		//display
	
	//----------------------------------------------------------------------  

}

void Tracer::writeImage(const std::string &fname) {
  CImg<unsigned char> output(res[0], res[1], 1, 3);
  for (unsigned int i=0; i<output.width(); i++) {
	for (unsigned int j=0; j<output.height(); j++) {
	  output(i, j, 0, 0) = std::min(1.0, std::max(0.0, im[j*output.width()+i][0]))*255.0;
	  output(i, j, 0, 1) = std::min(1.0, std::max(0.0, im[j*output.width()+i][1]))*255.0;
	  output(i, j, 0, 2) = std::min(1.0, std::max(0.0, im[j*output.width()+i][2]))*255.0;
	}
  }
  if (strstr(fname.c_str(), "png")) {
    output.save_png(fname.c_str());
	return;
  } else if (strstr(fname.c_str(), "jpg")) {
    output.save_jpeg(fname.c_str());
	return;
  }

#ifdef __APPLE__
  std::ofstream out(fname, std::ios::out | std::ios::binary);
#else
  std::ofstream out(fname.c_str(), std::ios_base::binary);
#endif
  out<<"P6"<<"\n"<<res[0]<<" "<<res[1]<<"\n"<<255<<"\n";
  Eigen::Vector3d *pixel = im;
  char val;
  for (unsigned int i=0; i<res[0]*res[1]; i++, pixel++) {
	val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[0])) * 255.0);
	out.write (&val, sizeof(unsigned char));
	val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[1])) * 255.0);
	out.write (&val, sizeof(unsigned char));
	val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[2])) * 255.0);
	out.write (&val, sizeof(unsigned char));
  }
  out.close();
}


int main(int argc, char *argv[]) {
  int c;
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
  while ((c = getopt(argc, argv, "a:s:d:cpjm:v")) != -1) {
	switch(c) {
	case 'a':
	  aperture = atof(optarg);
	  break;
	case 's':
	  samples = atoi(optarg);
	  break;
	case 'j':
	  stratified = true;
	  break;
	case 'c':
	  color = true;
	  break;
	case 'p':
	  phong = true;
	  break;
	case 'm':
	  mode = atoi(optarg);
	  break;
	case 'd':
	  maxraydepth = atoi(optarg);
	  break;
	case 'v':
	  verbose = true;
	  break;
	default:
	  abort();
	}
  }

  if (argc-optind != 2) {
	std::cout<<"usage: trace input.nff output.ppm"<<std::endl;
	for (int i=0; i<argc; i++) std::cout<<argv[i]<<std::endl;
	exit(0);
  }	

  switch(mode) {
  case 0:
    reflections = shadows = false;
    break;
  case 1:
    reflections = false;
    break;
  case 2:
    shadows = false;
    break;
  }

  Tracer tracer(argv[optind++]);
  tracer.aperture = aperture;
  tracer.samples = samples;
  tracer.stratified = stratified;
  tracer.color = color;
  tracer.phong = phong;
  tracer.reflections = reflections;
  tracer.shadows = shadows;
  tracer.maxraydepth = maxraydepth;
  tracer.verbose = verbose;
  tracer.createImage();
  std::cout<<"writing "<<argv[optind]<<std::endl;
  tracer.writeImage(argv[optind++]);
};
