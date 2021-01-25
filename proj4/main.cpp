#include "io.h"
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <sstream>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <cmath>

void laplacian(const std::vector<Tri> &tri, const std::vector<Eigen::Vector3d> &x, std::vector<Eigen::Vector3d> &y){
	std::vector<Eigen::Vector3d> lap;
	std::vector<double> m;
	lap.resize(x.size());
	m.resize(x.size());
	
	
	for(int j=0; j<x.size(); j++){
			m[j] = 0;
			lap[j] = {0.0, 0.0, 0.0};
		}
	for(int t = 0; t< tri.size(); t++){
		lap[tri[t][0]] += x[tri[t][1]] - x[tri[t][0]];
		lap[tri[t][0]] += x[tri[t][2]] - x[tri[t][0]];
		
		lap[tri[t][1]] += x[tri[t][0]] - x[tri[t][1]];
		lap[tri[t][1]] += x[tri[t][2]] - x[tri[t][1]];
		
		lap[tri[t][2]] += x[tri[t][0]] - x[tri[t][2]];
		lap[tri[t][2]] += x[tri[t][1]] - x[tri[t][2]];
			
		m[tri[t][0]] += 2;
		m[tri[t][1]] += 2;
		m[tri[t][2]] += 2;
	}
	for(int j = 0; j < x.size(); j++){
		y[j] = lap[j]/m[j];
	}
	
}


int main(int argc, char *argv[]){
	
	std::vector<Tri> tri;
	std::vector<Eigen::Vector3d> pts;
	std::vector<double> m;
	std::vector<Eigen::Vector3d> lap;
	
	readObjFile(argv[1], pts, tri);
	
	double stepSize = atoi(argv[3]);
	int numIterations = atof(argv[4]);
	
	for(int i= 0; i < numIterations; i++){
		//#IFDEF 0
	/*	for(int j=0; j<pts.size(); j++){
			m[j] = 0;
			lap[j] = 0;
		}
		for(int t = 0; t< (int) tri.size(); t++){
			lap[tri[t][0]] += pts[tri[t][1]] - pts[tri[t][0]];
			lap[tri[t][0]] += pts[tri[t][2]] - pts[tri[t][0]];
			
			lap[tri[t][1]] += pts[tri[t][0]] - pts[tri[t][1]];
			lap[tri[t][1]] += pts[tri[t][2]] - pts[tri[t][1]];
			
			lap[tri[t][2]] += pts[tri[t][0]] - pts[tri[t][2]];
			lap[tri[t][2]] += pts[tri[t][1]] - pts[tri[t][2]];
			
			m[tri[t][0]] += 2;
			m[tri[t][1]] += 2;
			m[tri[t][2]] += 2;
		}
		for(int j = 0; j < pts.size(); j++){
			pts[j] = pts[j] + atof(argv[4]) * lap[j]/m[j];
		}*/
//#else
	std::vector<Eigen::Vector3d> y;
	y.resize(pts.size());
	laplacian(tri, pts, y);
	for(int j = 0; j < pts.size(); j++){
		//pts[j] = pts[j] + atof(argv[4]) * y[j];
		pts[j] = pts[j] + stepSize * y[j];
		}
	}
	writeObjFile(argv[2], pts, tri);
}