/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	
	// TODO: Create normal distributions for y and theta
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	num_particles=120;
	
	
	for(unsigned int i=0;i<num_particles;i++)
	{
		Particle p;
		p.id=i;
		p.x=dist_x(gen);
		p.y=dist_y(gen);
		p.theta=dist_theta(gen);
		p.weight=1.0;
		particles.push_back(p);
		weights.push_back(1);
	

	}
	is_initialized=true;
	
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
		default_random_engine gen;
		//will be added as noise to the equations
	normal_distribution<double> dist_x(0, std_pos[0]); 
	
	// TODO: Create normal distributions for y and theta
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	for (int i=0;i<particles.size();i++)
	{
		Particle &p=particles[i];
		
		
		if((yaw_rate<0.0001 && yaw_rate>-0.0001))
		{
		p.x+=(velocity)*delta_t*(cos(p.theta))+dist_x(gen);
		p.y+=(velocity)*delta_t*(sin(p.theta))+dist_y(gen);
		p.theta+=dist_theta(gen);
		
		}
		
		else
		{
		
		p.x+=(velocity/yaw_rate)*(sin(p.theta+yaw_rate*delta_t)-sin(p.theta))+dist_x(gen);
		p.y+=(velocity/yaw_rate)*(cos(p.theta)-cos(p.theta+yaw_rate*delta_t))+dist_y(gen);
		p.theta+=yaw_rate*delta_t+dist_theta(gen);
		}
		
	}
	
	
	
}

void data_association(std::vector<LandmarkObs> &predicted, std::vector<LandmarkObs>& observations,std::vector<LandmarkObs>& particle_lms) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	
	
	double min_dist;
	LandmarkObs best_lm;
	for (int i=0;i<observations.size();i++)
	{
		
		min_dist=sqrt((observations[i].x-particle_lms[0].x)*(observations[i].x-particle_lms[0].x)+(observations[i].y-particle_lms[0].y)*(observations[i].y-particle_lms[0].y));
		
		for(int j=0; j< particle_lms.size(); j++){

			
			double dist = sqrt((particle_lms[j].x - observations[i].x)*(particle_lms[j].x - observations[i].x) + (particle_lms[j].y - observations[i].y)*(particle_lms[j].y - observations[i].y));
			
			
			if(dist < min_dist ){
				min_dist = dist;
				best_lm = particle_lms[j]; //closest landmark
			}
		}
		
		 
			//cout<<"Closest landmark is"<<predicted_landmark_ind;
		predicted.push_back(best_lm);
				
	}
	
	
	
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> &observations, Map &map_landmarks) {
			
		double exp1 = 1.0/(2 * M_PI * std_landmark[0] * std_landmark[1]);
		double power_stdx=std_landmark[0]*std_landmark[0];
		double power_stdy=std_landmark[1]*std_landmark[1];
		
	for (int i = 0; i < num_particles; i++) {
		
		vector<LandmarkObs> trans_lms ;
		
		
		// For each landmark transform the landmark in the particle system and store in a transformed landmark list
		for (int j = 0; j< map_landmarks.landmark_list.size(); j++)
		{
			Map::single_landmark_s each_landmark = map_landmarks.landmark_list[j];
			
				int id = each_landmark.id_i;
				double range_x=each_landmark.x_f - particles[i].x;
				double range_y=each_landmark.y_f - particles[i].y;
				double x_l = range_x * cos(particles[i].theta) + range_y * sin(particles[i].theta); // cos(theta-pi/2)=sin(theta) and sin(theta-pi/2)=-cos(theta)
				double y_l = -range_x * sin(particles[i].theta) + range_y * cos(particles[i].theta);

				
				
				//And add to the a  new landmark list
				trans_lms.push_back(LandmarkObs{id,x_l,y_l});

		}
			
			// Now associate observations with the transformed landmark list  
			std::vector<LandmarkObs> predicted ;
			data_association(predicted, observations, trans_lms);
			
		
		//Calculate multi-variate guassian on each observation and its closest landmark 
		
		
		double weight_new = 0;

		
		for(int k=0; k< observations.size();k++){
		//Compute predicted landmark measurement for the particle

		double xDiff = observations[k].x - predicted[k].x;
		double yDiff = observations[k].y - predicted[k].y;

		double power_xDiff=xDiff*xDiff;
		double power_yDiff=yDiff*yDiff;
		double exp2 = exp(-0.5 * (power_xDiff/power_stdx + power_yDiff/power_stdy));
		double weight_each = exp1*exp2;
		
		//multiply density for all predicted measurements
		//for first measurement set product as the first weight
		if(weight_new != 0){
			weight_new = weight_new *weight_each;
		}else{
			weight_new=weight_each;
		}
		
	}
	//Update particle weight
	particles[i].weight = weight_new;
	weights[i] = weight_new;
	
}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	vector<Particle> Particles_re;
	// Take a discrete distribution with pmf equal to weights
    discrete_distribution<> resample_weights(weights.begin(), weights.end());
    // initialise new particle array
    
    // resample particles
    for (int i = 0; i < num_particles; ++i)
        Particles_re.push_back(particles[resample_weights(gen)]);

	particles = Particles_re;

}

Particle ParticleFilter::SetAssociations(Particle& particle,  std::vector<int>& associations, 
                                      std::vector<double>& sense_x,  std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
	
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
