
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>
#include <Eigen/Dense>
#include <cassert>
#include "kdTree.h"

using namespace std;
//using namespace cimg_library;


inline double sqr(double x) {return x*x;}
inline double cube(double x) {return x*x*x;}

//holds all of the boids properties from the file
struct boidProps{
	double size; //size of boid
	double nei_rad; //distance boid can see
	double mass_of_boid; //mass of boid
	double collision; //weights for forces
	double centering; //weights for forces
	double velocity; //weights for forces
	double hunger; //for extra credit
	double damping; //a constant
	double dt;  //timestep size
	double len_of_anim;  //length of animation
	int num_nei; //number of neighbors the boid can process
};

struct Boid{
	Eigen::Vector3d pos, vel, force;
	void move(const boidProps &boidProps);
};


void Boid::move(const boidProps &boidProps){
	vel += boidProps.dt * force/boidProps.mass_of_boid;
	vel *= boidProps.damping;
	pos += boidProps.dt * vel;
	
	if(pos[0] < -0.5){
		pos[0] = -0.5;
		vel[0] = fabs(vel[0]); //fabs is floating point absolute value
	}
	if(pos[1] < -0.25){
		pos[1] = -0.25;
		vel[1] = fabs(vel[1]); //fabs is floating point absolute value
	}
	if(pos[2] < -0.125){
		pos[2] = -0.125;
		vel[2] = fabs(vel[2]); //fabs is floating point absolute value
	}
	if(pos[0] > 0.5){
		pos[0] = 0.5;
		vel[0] = -fabs(vel[0]); //fabs is floating point absolute value
	}
	if(pos[1] > 0.25){
		pos[1] = 0.25;
		vel[1] = -fabs(vel[1]); //fabs is floating point absolute value
	}
	if(pos[2] > 0.125){
		pos[2] = 0.125;
		vel[2] = -fabs(vel[2]); //fabs is floating point absolute value
	}
	
}


void timestep(const boidProps &boidProps, std::vector<Boid> &boids){
	std::vector<Eigen::Vector3d> points;
	points.resize(boids.size());
	for(unsigned int i=0; i<boids.size(); i++){
		points[i] = boids[i].pos;
	}


	KDTree tree(points);
	for(unsigned int i=0; i<boids.size(); i++){
		std::vector<int> neighbors;
		Boid &boid = boids[i];
		tree.neighbors(points, boid.pos, boidProps.num_nei, boidProps.nei_rad, neighbors);
		boid.force<<0.0,0.0,0.0;
		//average position and average velocity
		Eigen::Vector3d avgpos, avgvel;
		avgpos<<0.0,0.0,0.0;
		avgvel<<0.0,0.0,0.0;
		//counter for neighbors
		int c = 0;
		for(unsigned int j=0; j<neighbors.size(); j++){
			int n = neighbors[j];
			if(n = i){
				continue;
			}
			c++;
			Eigen::Vector3d force = boid.pos - boids[n].pos;
			double m = force.norm();
			if(m > 0.0001){
				boid.force += boidProps.collision * force/cube(m);
			}
			avgpos += boids[n].pos;
			avgvel += boids[n].vel;
		}
		//assert(c == neighbors.size()-1);
		avgpos /= c;
		avgvel /= c;
		Eigen::Vector3d center = avgpos - boid.pos;
		double m = center.norm();
		if(m > 0.0001){
			boid.force += boidProps.centering * (avgpos - boid.pos)/m;
		}
	}
	for(int i = 0; i < boids.size(); i++){
		boids[i].move(boidProps);
	}
	
}


int main(int argc, char *argv[]) {
	if(argc > 1){
	std::ifstream in(argv[1], std::ios_base::in);
	std::string line;
	char ch;
	boidProps props;

	int num_of_boids;
		
	vector<Boid> boids;
 //parser
	in>>props.size>>
	props.nei_rad>>
	props.num_nei>>
	props.mass_of_boid>>
	props.collision>>
	props.centering>>
	props.velocity>>
	props.hunger>>
	props.damping>>
	props.dt>>
	props.len_of_anim;
	/*
	cout<<props.size<<" "<<
	ch<<props.nei_rad<<" "<<
	ch<<props.num_nei<<" "<<
	ch<<props.mass_of_boid<<" "<<
	ch<<props.collision<<" "<<
	ch<<props.centering<<" "<<
	ch<<props.velocity<<" "<<
	ch<<props.hunger<<" "<<
	ch<<props.damping<<" "<<
	ch<<props.dt<<" "<<
	ch<<props.len_of_anim<<endl;
	*/
	getline(in, line);
	in>>num_of_boids;
	
	cout<<num_of_boids<<endl;

	for(int i = 0; i < num_of_boids; i++){ 
		Boid boid;
		in>>ch;
		getline(in, line, ',');
		std::stringstream ss(line);
		ss>>boid.pos[0];
		getline(in, line, ',');
		ss= std::stringstream(line);
		ss>>boid.pos[1];
		getline(in, line, ']');
		ss= std::stringstream(line);
		ss>>boid.pos[2];
		
		getline(in, line, '[');
		
		getline(in, line, ',');
		ss= std::stringstream(line);
		ss>>boid.vel[0];
		getline(in, line, ',');
		ss= std::stringstream(line);
		ss>>boid.vel[1];
		getline(in, line, ']');
		ss= std::stringstream(line);
		ss>>boid.vel[2];
		
		boids.push_back(boid);
		
		//cout<<boid.pos[0]<<" "<<boid.pos[1]<<" "<<boid.pos[2]<<" "<<boid.vel[0]<<" "<<boid.vel[1]<<" " <<boid.vel[2]<< endl;
	}
//end parser

//making the outfile
	std::ofstream test(argv[2]);
	
	int frames = (int) ceil(props.len_of_anim/props.dt);
	test<<frames<<endl;
	
	for(double t = 0; t < props.len_of_anim; t+= props.dt){
		timestep(props, boids);
		test<<boids.size()<<endl;
		for(int i = 0; i < boids.size(); i++){
		test<<"["<<boids[i].pos[0]<<","<<boids[i].pos[1]<<","<<boids[i].pos[2]
		<<"] ["<<
		boids[i].vel[0]<<","<<boids[i].vel[1]<<","<<boids[i].vel[2]<<"]"<<endl;
		}
		test<<0<<endl;
	} 
	   
	   
	}
	return 0;
}; 