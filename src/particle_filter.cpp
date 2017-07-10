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
#include <initializer_list>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  cout << "debug: Enter init" << endl;
  default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  num_particles = 1000;
  for (int i=0; i<num_particles; ++i){
    Particle particle;
    particle = {i, dist_x(gen), dist_y(gen), dist_theta(gen), 1};
    particles.push_back(particle);
    weights.push_back(1);
  }
  is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    cout << "debug: Enter prediction" << endl;
    default_random_engine gen;

    for (int i=0; i < particles.size(); ++i){
        normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
        normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
        normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);


        double x = dist_x(gen);
        double y = dist_y(gen);
        double theta = dist_theta(gen);
        cout << "debug: x " << x << endl;
        cout << "debug: y " << y << endl;
        cout << "debug: theta " << theta << endl;
        if (fabs(yaw_rate) >= 0.001){
            particles[i].x = x + velocity/yaw_rate *(sin(theta + yaw_rate*delta_t) - sin(theta));
            particles[i].y = y + velocity/yaw_rate *(- cos(theta + yaw_rate*delta_t) + cos(theta));
            particles[i].theta = theta + yaw_rate*delta_t;
        }
        else{
            particles[i].x = x + velocity * delta_t * cos(theta);
            particles[i].y = y + velocity * delta_t * sin(theta);
            particles[i].theta = theta;
        }
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    for (unsigned int i=0; i<observations.size(); ++i){
        double min_dist = numeric_limits<double>::max();
        double obs_x = observations[i].x;
        double obs_y = observations[i].y;
        for (unsigned int j=0; j<predicted.size(); ++j){
            double pred_x = predicted[j].x;
            double pred_y = predicted[j].y;
            double dist = pow(pred_x-obs_x, 2)+pow(abs(pred_y-obs_y), 2);
            if (dist < min_dist) {
                observations[i].id = predicted[j].id;
                min_dist = dist;
            }
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    cout << "debug: Enter update weights" << endl;
    for (unsigned int i=0; i<particles.size(); ++i){
        double p_x = particles[i].x;
        double p_y = particles[i].y;
        double p_theta = particles[i].theta;

        // transform observations from VEH to MAP coordinates
        std::vector<LandmarkObs> obs_transform(observations.size());
        for (unsigned int j=0; j<observations.size();j++){
            double obs_x = observations[j].x;
            double obs_y = observations[j].y;
            obs_transform[j].x = p_x + cos(p_theta)*obs_x - sin(p_theta)*obs_y;
            obs_transform[j].y = p_y + sin(p_theta)*obs_x + cos(p_theta)*obs_y;
        }

        // get landmarks within sensor range
        std::vector<LandmarkObs> predicted;
        for (unsigned int j=0; j<map_landmarks.landmark_list.size(); ++j){
            double lm_x = map_landmarks.landmark_list[j].x_f;
            double lm_y = map_landmarks.landmark_list[j].y_f;
            int lm_id = map_landmarks.landmark_list[j].id_i;
            if (pow(lm_x-p_x,2) + pow(lm_y-p_y,2) < sensor_range){
                LandmarkObs lm = {lm_id, lm_x, lm_y};
                predicted.push_back(lm);
            }
        }
        // Assign closest landmarks to observation via id
        dataAssociation(predicted, obs_transform);

        vector<int> associations;
        vector<double> sense_x;
        vector<double> sense_y;
        for (unsigned int j=0; j<obs_transform.size(); ++j){
            associations.push_back(obs_transform[j].id);
            sense_x.push_back(obs_transform[j].x);
            sense_y.push_back(obs_transform[j].y);
        }
        particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);

        // Update weights
        particles[i].weight = 1./2/M_PI/std_landmark[0]/std_landmark[1];
        for (unsigned int j=0; j<obs_transform.size(); ++j){
            for (unsigned int k=0; k<predicted.size(); ++k){
                if (obs_transform[j].id == predicted[k].id){
                    double obs_x = obs_transform[j].x;
                    double obs_y = obs_transform[j].y;
                    double mu_x = predicted[k].x;
                    double mu_y = predicted[k].y;
                    particles[i].weight *= exp(-0.5*(pow(obs_x-mu_x,2)*pow(std_landmark[0],2) +
                                  pow(obs_y-mu_y,2)*pow(std_landmark[1],2)));
                }
            }

        }
    }
    for (unsigned int i=0; i<particles.size(); ++i){
        weights[i] = particles[i].weight;
        //cout << "debug: " << weights[i] << endl;
    }
    //for (int i=0; i<num_particles;i+=100){
    //    cout << "debug: weights" << weights[i] << endl;
    //}

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
    // NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    cout << "debug: Enter resample" << endl;

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> d(weights.begin(), weights.end());
    //particles.clear();
    vector<Particle> particles_resampled;
    for (unsigned int i=0; i<num_particles; ++i){
        Particle p;
        int rand_int = d(gen);
        p.id = i;
        p.x = particles[rand_int].x;
        p.y = particles[rand_int].y;
        p.theta = particles[rand_int].theta;
        p.weight = particles[rand_int].weight;
        particles_resampled.push_back(p);
    }
    particles = particles_resampled;


}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

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
