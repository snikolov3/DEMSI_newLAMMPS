#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "random_mars.h"
#include "lmp_data.h"

void gen_nonoverlapping(const std::vector<double> bin_edges, const std::vector<double> bin_values,
    int interpolate, int seed, int maxattempts, double vf_target, Lmp_data *lmp);

template <typename T>
T lexical_cast(const std::string& str)
{
    T out;
    std::istringstream iss;
    iss.str(str);
    iss >> out;
    return out;
}

bool is_digits(const std::string &str)
{
    return str.find_first_not_of("0123456789.") == std::string::npos;
}

int main(int argc, char** argv){
  std::ifstream infile;
  std::vector <Particle> particles;

  int valflag;
  std::vector <double> bin_edges;
  std::vector <double> bin_values;
  std::vector <std::string> words;
  std::string str;

  int ndim;
  int seed;
  double xlo, xhi, ylo, yhi;
  int xperiodic, yperiodic, zperiodic;
  double vf_target;
  int interpolate;
  int maxattempts = 100;

  //Parse args
  if (argc != 14){
    std::cout<<"Usage: "<< argv[0] <<" <name of file containing particle size distribution> <dimensionality (2 or 3)> "
        "<xlo> <xhi> <ylo> <yhi> "
        "<xperiodic> <yperiodic> "
        "<interpolate flag (0 or 1)> <target volume fraction> "
        "<max attempts per particle> <random seed> <output file name>"<<std::endl;
    return 0;
  }
  ndim = atoi(argv[2]);
  xlo = atof(argv[3]);
  xhi = atof(argv[4]);
  ylo = atof(argv[5]);
  yhi = atof(argv[6]);
  //zlo = atof(argv[7]);
  // zhi = atof(argv[8]);
  xperiodic = atoi(argv[7]);
  yperiodic = atoi(argv[8]);
  zperiodic = 1;
  //zperiodic = atoi(argv[11]);
  interpolate = atoi(argv[9]);
  vf_target = atof(argv[10]);
  maxattempts = atoi(argv[11]);
  seed = atoi(argv[12]);

  //Parse PSD file, populate histogram
  valflag = 0;
  infile.open(argv[1], std::ios::in);
  if (infile.is_open()){
    while (std::getline(infile, str)){
      if (!is_digits(str) || str==""){
        valflag = 1;
        continue;
      }

      if (valflag) bin_values.push_back(lexical_cast<double>(str));
      else bin_edges.push_back(lexical_cast<double>(str));
    }
  }
  else{
    std::cout <<"Failed to open file "<<argv[1]<<std::endl;
    return -1;
  }
  //Add error check that binvalues/edges have the correct dimensions

  double maxdiam = 0;
  for (int i = 0; i < bin_edges[i]; i++){
    if (bin_edges[i] > maxdiam) maxdiam = bin_edges[i];
  }

  //Create instance of Lmp_data
  Lmp_data lmp(0, ndim, xperiodic, yperiodic, zperiodic, xlo, xhi, ylo, yhi, -maxdiam, maxdiam);

  //Run particle placement
  gen_nonoverlapping(bin_edges, bin_values, interpolate, seed, maxattempts, vf_target, &lmp);
  //Output to lammps data file format
  lmp.write_data(argv[13]);


}



void gen_nonoverlapping(const std::vector<double> bin_edges, std::vector<double> bin_values,
    int interpolate, int seed, int maxattempts, double vf_target, Lmp_data *lmp)
{
  const double PI  =3.141592653589793238463;
  int nattempts;
  double draw, rad;
  double x, y, z;
  std::vector <double> radvals;

  double vol = 0;
  double vf;
  int ibin;

  int nbins = (int)bin_values.size();

  double cdf[nbins+1];
  double maxrad = 0;

  //Compute (normalized) c.d.f
  double norm = 0;
  for (int i = 0; i < nbins; i++){
    norm += (bin_edges[i+1]-bin_edges[i])*bin_values[i];
  }
  for (int i = 0; i < nbins+1; i++){
    cdf[i] = 0;
    for (int j = 1; j <= i; j++){
      cdf[i] += (bin_edges[j]-bin_edges[j-1])*bin_values[j-1]/norm;
    }
  }
  //for (int i = 0; i < nbins+1; i++) cdf[i] /= cdfsum;

  RanMars* randgen = new RanMars(seed);
  lmp->ntot = 0;

  //Draw particle radii from distribution
  while (1){
    draw = randgen->uniform();
    for (ibin = 0; ibin < nbins; ibin++){
      if (draw < cdf[ibin]) break;
    }
    if (interpolate){
      double frac = (draw - cdf[ibin-1])/(cdf[ibin]-cdf[ibin-1]);
      rad = bin_edges[ibin-1] + frac*(bin_edges[ibin]-bin_edges[ibin-1]);
    }
    else rad = 0.5*(bin_edges[ibin]+bin_edges[ibin-1]);
    rad *= 0.5; //assuming file contains diameters
    radvals.push_back(rad);

    if (lmp->ndim == 2){
      vol += PI*rad*rad;
      vf = vol/(lmp->lx*lmp->ly);
    }
    else{
      vol += 4.0/3.0*PI*rad*rad*rad;
      vf = vol/(lmp->lx*lmp->ly*lmp->lz);
    }
    if (vf > vf_target) break;
  }

  std::cout << "Attempting to place "<< radvals.size() << " particles." << std::endl;

  //Sort from largest to smallest
  std::sort(radvals.rbegin(), radvals.rend());
  maxrad = radvals[0];
  if (lmp->ndim == 2){
    lmp->zlo = -1.1*maxrad;
    lmp->zhi = 1.1*maxrad;
    lmp->lz = lmp->zhi - lmp->zlo;
  }

  lmp->particles = new Particle[(int)radvals.size()];

  //Since maxrad was previously 0, bins were not set up, set them up here
  lmp->setup_bins(maxrad);

  //Place particles
  int nfailed = 0;
  int id = 1;
  vol = 0;
  for (int i = 0; i < (int)radvals.size(); i++){
    rad = radvals[i];
    nattempts = 0;
    if (i%10 == 0) std::cout << "Attempting to place particle " << lmp->ntot <<
        " with radius "<<rad << std::endl;
    while(1){
      nattempts++;
      if (nattempts > maxattempts){
        std::cout << "Failed to place particle with radius "<<rad << std::endl;
        nfailed++;
        break;
      }
      x = randgen->uniform()*lmp->lx;
      y = randgen->uniform()*lmp->ly;
      z = 0;
      if (lmp->ndim == 3)
        z = randgen->uniform()*lmp->lz;

      if (!lmp->xperiodic){
        if (x > lmp->lx - rad || x < rad) {
          continue;
        }
      }
      if (!lmp->yperiodic){
        if (y > lmp->ly - rad || y < rad) {
          continue;
        }
      }
      if (lmp->ndim ==3 && !lmp->zperiodic){
        if (z > lmp->lz - rad || z < rad) {
          continue;
        }
      }
      if (!lmp->check_overlap(x, y, z, rad)){
        lmp->particles[lmp->ntot].id = id++;
        lmp->particles[lmp->ntot].x = x;
        lmp->particles[lmp->ntot].y = y;
        lmp->particles[lmp->ntot].z = z;
        lmp->particles[lmp->ntot].rad = rad;
        lmp->particles[lmp->ntot].diam = 2.0*rad;
        lmp->particles[lmp->ntot].type = 1;
        lmp->particles[lmp->ntot].density = 1.0;
        lmp->add_particle_to_bins(lmp->particles[lmp->ntot], 2.0*maxrad);
        lmp->ntot++;
        if (lmp->ndim == 2){
          vol += PI*rad*rad;
          vf = vol/(lmp->lx*lmp->ly);
        }
        else{
          vol += 4.0/3.0*PI*rad*rad*rad;
          vf = vol/(lmp->lx*lmp->ly*lmp->lz);
        }
        if (i%10 == 0) std::cout << "Particle placed successfully, volume fraction now "<< vf <<
            " target is "<< vf_target << std::endl;
        break;
      }     
    }
    if (vf > vf_target) break;
  }
  std::cout << "Failed to place a total of "<<nfailed<<" particles"<<std::endl;
}
