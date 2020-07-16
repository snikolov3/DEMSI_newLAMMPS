#ifndef LMP_DATA_H
#define LMP_DATA_H

#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>

#define MAX
struct Particle{
	int id, molid, type;
	double diam, rad, rad2, density;
	double x, y, z;
	int ix, iy, iz;
};

struct Bininfo{
  int n;
  std::vector < Particle > list;
};

class Lmp_data{
public:
  Lmp_data(int, int,
        int, int, int,
        double, double, double, double, double, double);
  virtual ~Lmp_data();

  int ndim;

  int ntot;
  int xperiodic, yperiodic, zperiodic;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double lx, ly, lz;

  double minrad, maxrad;
  Particle *particles;

  int shiftflag;
  int nbinx, nbiny, nbinz;
  double binsizex, binsizey, binsizez;
  Bininfo ***bins;

  void setup_bins(double);
  void bin_all_particles(int);
  void add_particle_to_bins(Particle, double);
  int point_in_sphere(const double, const double, const double, const Particle);
  int check_overlap(double, double, double, double);

  void write_data(const char*);
private:
  double bindiag;
};

#endif


