/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;

// declare a random engine to be used across multiple and various method calls
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  std::normal_distribution<double> noise_x(0,std[0]);
  std::normal_distribution<double> noise_y(0,std[1]);
  std::normal_distribution<double> noise_theta(0,std[2]);

  /*Random seed*/
  std::random_device rd;
  std::default_random_engine generator( rd() );

  for (int i = 0; i < num_particles; ++i)
  {
    /* code */
    Particle p;

    p.id = i;
    p.x = x + noise_x(generator);
    p.y = y + noise_y(generator);
    p.theta = theta + noise_theta(generator);
    p.weight = 1;

    particles.push_back(p);
  }

  is_initialized = true;
  

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::normal_distribution<double> noise_x(0,std_pos[0]);
  std::normal_distribution<double> noise_y(0,std_pos[1]);
  std::normal_distribution<double> noise_theta(0,std_pos[2]);

  /*Random seed*/
  std::random_device rd;
  std::default_random_engine generator( rd() );


   for (int i = 0; i < num_particles; ++i)
   {
    if(yaw_rate != 0.0 ){
      particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
      particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      particles[i].theta += yaw_rate * delta_t;
    }
    else{
      particles[i].x += velocity * cos(particles[i].theta) * delta_t; 
      particles[i].y += velocity * sin(particles[i].theta) * delta_t; 
    }

    /*Adding noise*/
    particles[i].x += noise_x(generator);
    particles[i].y += noise_y(generator);
    particles[i].theta += noise_theta(generator);

  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  for (long unsigned int i = 0; i < observations.size(); ++i)
  {
    double minimum_dist = std::numeric_limits<double>::max();
    int obs_id = -1;

    for (long unsigned int j = 0; j < predicted.size(); ++j)
    {
      double actual_dist = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);

      if(actual_dist < minimum_dist){
        minimum_dist = actual_dist;
        obs_id = predicted[j].id;
      }

    }

    observations[i].id = obs_id;

  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  for (long unsigned int i = 0; i < particles.size(); ++i)
  {
    double particles_x = particles[i].x;
    double particles_y = particles[i].y;
    double particles_theta = particles[i].theta;


    /*transforming sensor observations to map observations*/
    vector<LandmarkObs> map_observations;
    for (long unsigned int j = 0; j < observations.size(); ++j)
    {
      LandmarkObs temp_obs;
      temp_obs.id = observations[j].id;     
      temp_obs.x = particles_x + (cos(particles_theta) * observations[j].x) - (sin(particles_theta) * observations[j].y);   
      temp_obs.y = particles_y + (sin(particles_theta) * observations[j].x) + (cos(particles_theta) * observations[j].y); 
      map_observations.push_back(temp_obs);
    }

    /*setting predictions*/
    vector<LandmarkObs> predictions;
    for (long unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j)
    {
      double landmark_x = map_landmarks.landmark_list[j].x_f;
      double landmark_y = map_landmarks.landmark_list[j].y_f;
      int landmark_id = map_landmarks.landmark_list[j].id_i;

      if (fabs(landmark_x - particles_x) <= sensor_range && fabs(landmark_y - particles_y) <= sensor_range)/*fast algorithm*/
      /*double dist2point = dist(landmark_x, landmark_y, particles_x, particles_y);
      if (dist2point <= sensor_range)*/
      {
        predictions.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
      }

    }

    /*Find the predicted measurement that is closest to each observed measurement and assign the observed measurement to this particular landmark.*/
    dataAssociation(predictions,map_observations);

    /*Reseting weight to 1*/
    particles[i].weight = 1.0;

    for (long unsigned int j = 0; j < map_observations.size(); ++j)
    {

      double obs_x = map_observations[j].x;
      double obs_y = map_observations[j].y;
      double pre_x = 0;
      double pre_y = 0;

      /*finding the predicted point that match the observations for multivariate Gaussian calculation*/
      for (long unsigned int k = 0; k < predictions.size(); ++k)
      {
        if (map_observations[j].id == predictions[k].id){
          pre_x = predictions[k].x;
          pre_y = predictions[k].y;
          break;
        }
      }
      
      /*Weight calculation using multivariate Gaussian*/
      double s_l_x = std_landmark[0];
      double s_l_y = std_landmark[1];
      double obs_weight = ( 1/(2*M_PI*s_l_x*s_l_y)) * exp( (-0.5)*( pow(pre_x-obs_x,2)/(pow(s_l_x, 2)) + (pow(pre_y-obs_y,2)/(pow(s_l_y, 2))) ) );

      particles[i].weight *= obs_weight;
    }

  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  /*Random seed*/
  std::random_device rd;
  std::default_random_engine generator( rd() );

  weights.erase (weights.begin(),weights.end());

  for (long unsigned int i = 0; i < particles.size(); ++i)
  {
    weights.push_back(particles[i].weight);
  }

  /*create new particle vector*/
  vector<Particle> new_particles;

  /*get the max value of the weight in the particles vector*/
  double wmax = *std::max_element(weights.begin(),weights.end());
  /*create a uniform real distribution from 0.0 to wmax*/
  std::uniform_real_distribution<double> uniform_real_dist_weights(0.0, wmax);

  /*create a uniform distribution from 0 to (number of particles -1)*/
  std::uniform_int_distribution<int> uni_int_dist_index(0, num_particles-1);
  int index = uni_int_dist_index(generator);

  double beta = 0.0;

  for (int i = 0; i < num_particles; ++i)
  {
    beta += 2.0 * uniform_real_dist_weights(generator);
    while(beta > weights[index]){
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }

  particles = new_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}