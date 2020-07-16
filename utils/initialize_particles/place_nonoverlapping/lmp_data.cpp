#include "lmp_data.h"

Lmp_data::Lmp_data(int bintype, int ndim_, int xperiodic_, int yperiodic_, int zperiodic_,
    double xlo_, double xhi_, double ylo_, double yhi_, double zlo_, double zhi_):
    ndim(ndim_), xperiodic(xperiodic_), yperiodic(yperiodic_), zperiodic(zperiodic_),
    xlo(xlo_), xhi(xhi_), ylo(ylo_), yhi(yhi_), zlo(zlo_), zhi(zhi_)
{

  lx = xhi - xlo;
  ly = yhi - ylo;
  lz = zhi - zlo;

  maxrad = -1;
  minrad = 1e20;
  shiftflag = 1;

  if (maxrad > 0){
    setup_bins(maxrad);
    //Bintype = 0 - set up bins, don't fill them
    //Bintype = 1 - radius for bin coverage is radius of each particle
    //Bintype = 2 - radius for bin coverage is 2*std::max. radius (if particle contacts are later searched)
    printf("bintype is %d\n",bintype);
    if (bintype > 0){
      printf("Binning %d particles...",ntot);
      bin_all_particles(bintype);
      printf("done.\n");
    }
  }
}

Lmp_data::~Lmp_data(){
  if (shiftflag){
    for (int i = 0; i < ntot; i++){
      particles[i].x += xlo;
      particles[i].y += ylo;
      if (ndim == 3) particles[i].z += zlo;
    }
  }
  for (int i = 0; i < nbinx; i++){
    for (int j = 0; j < nbiny; j++){
      delete [] bins[i][j];
    }
    delete [] bins[i];
  }
  delete [] bins;
}

void Lmp_data::setup_bins(double maxrad){
  if (ndim == 2){
    binsizex = std::min(lx,std::max(maxrad,lx/1000.0));
    binsizey = std::min(ly,std::max(maxrad,ly/1000.0));
    binsizez = lz;
  }
  else{
    binsizex = std::min(lx,std::max(maxrad,lx/100.0));
    binsizey = std::min(ly,std::max(maxrad,ly/100.0));
    binsizez = std::min(lz,std::max(maxrad,lz/100.0));
  }

  nbinx = static_cast<int>(lx/binsizex);
  nbiny = static_cast<int>(ly/binsizey);
  nbinz = static_cast<int>(lz/binsizez);

  binsizex = lx/nbinx;
  binsizey = ly/nbiny;
  binsizez = lz/nbinz;

  //In case particles are exactly at xhi, yhi, zhi
  nbinx+=1;
  nbiny+=1;
  if (ndim == 3) nbinz+=1;

  bins = new Bininfo**[nbinx];
  for (int i = 0; i < nbinx; i++){
    bins[i] = new Bininfo*[nbiny];
    for (int j = 0; j < nbiny; j++){
      bins[i][j] = new Bininfo[nbinz];
      for (int k = 0; k < nbinz; k++)
        bins[i][j][k].n = 0;
    }
  }
}

void Lmp_data::bin_all_particles(int bintype){
  double binrad;
  //Put particles in lists for each bin
  for (int i = 0; i < ntot; i++){
    if (bintype == 1) binrad = particles[i].rad;
    else if (bintype == 2) binrad = 2.0*maxrad;
    add_particle_to_bins(particles[i], binrad);
  }
}


int Lmp_data::point_in_sphere(const double x, const double y, const double z, const Particle p){
  double dx, dy, dz;
  dx = p.x - x;
  dy = p.y - y;
  dz = p.z - z;
  return (dx*dx + dy*dy + dz*dz < p.rad*p.rad);
}

void Lmp_data::add_particle_to_bins(Particle p, double binrad){
  int ixlo, ixhi, iylo, iyhi, izlo, izhi;
  int iix, iiy, iiz, ix, iy, iz;
  double xs, ys, zs, xc, yc, zc;

  double diagsq = (binsizex*binsizex + binsizey*binsizey + binsizez*binsizez);
  bindiag = sqrt(diagsq);

  Particle s;

  ixlo = static_cast<int>((p.x - binrad)/binsizex)-1;
  if (p.x < binrad) ixlo -= 1; //Because e.g. -0.5 becomes 0
  ixhi = static_cast<int>((p.x + binrad)/binsizex)+1;
  if (ixlo < 0 && !xperiodic) ixlo = 0;
  if (ixhi >= nbinx && !xperiodic) ixhi = nbinx-1;

  iylo = static_cast<int>((p.y - binrad)/binsizey)-1;
  if (p.y < binrad) iylo -= 1;
  iyhi = static_cast<int>((p.y + binrad)/binsizey)+1;
  if (iylo < 0 && !yperiodic) iylo = 0;
  if (iyhi >= nbiny && !yperiodic) iyhi = nbiny-1;

  izlo = static_cast<int>((p.z - binrad)/binsizez)-1;
  if (p.z < binrad) izlo -= 1;
  izhi = static_cast<int>((p.z + binrad)/binsizez)+1;
  if (izlo < 0 && !zperiodic) izlo = 0;
  if (izhi >= nbinz && !zperiodic) izhi = nbinz-1;

  for (ix = ixlo; ix <= ixhi; ix++){
    iix = ix;
    xs = p.x;
    if (ix < 0 && xperiodic){
      iix = std::max(nbinx + ix,0);
      xs = p.x + lx;
    }
    if (ix >= nbinx && xperiodic){
      iix = std::min(ix - nbinx,nbinx-1);
      xs = p.x - lx;
    }
    for (iy = iylo; iy <= iyhi; iy++){
      iiy = iy;
      ys = p.y;
      if (iy < 0 && yperiodic){
        iiy = std::max(nbiny + iy,0);
        ys = p.y + ly;
      }
      if (iy >= nbiny && yperiodic){
        iiy = std::min(iy - nbiny, nbiny-1);
        ys = p.y - ly;
      }
      for (iz = izlo; iz <= izhi; iz++){
        iiz = iz;
        zs = p.z;
        if (iz < 0 && zperiodic){
          iiz = std::max(nbinz + iz,0);
          zs = p.z + lz;
        }
        if (iz >= nbinz && zperiodic){
          iiz = std::min(iz - nbinz,nbinz-1);
          zs = p.z - lz;
        }
        xc = (iix + 0.5)*binsizex;
        yc = (iiy + 0.5)*binsizey;
        zc = (iiz + 0.5)*binsizez;

        if ((xs-xc)*(xs-xc) + (ys-yc)*(ys-yc) + (zs-zc)*(zs-zc) < (binrad + bindiag)*(binrad + bindiag)){
          s.x = xs;
          s.y = ys;
          s.z = zs;
          s.type = p.type;
          s.rad = p.rad;
          s.rad2 = p.rad*p.rad;
          s.id = p.id;
          bins[iix][iiy][iiz].list.push_back(s);
          bins[iix][iiy][iiz].n++;
        }
      }
    }
  }

}

int Lmp_data::check_overlap(double x, double y, double z, double rad){
  double rx, ry, rz, r2;
  int ix = static_cast<int>(x/binsizex);
  int iy = static_cast<int>(y/binsizey);
  int iz = static_cast<int>(z/binsizez);
  Bininfo bin;

  bin = bins[ix][iy][iz];
  for (int i = 0; i < bin.n; i++){
    rx = x - bin.list[i].x;
    if (xperiodic) rx = rx - static_cast<int>(2.0*rx/lx)*lx;
    ry = y - bin.list[i].y;
    if (yperiodic) ry = ry - static_cast<int>(2.0*ry/ly)*ly;
    rz = z - bin.list[i].z;
    if (zperiodic) rz = rz - static_cast<int>(2.0*rz/lz)*lz;
    r2 = sqrt(rx*rx + ry*ry + rz*rz);
    if (r2 < rad + bin.list[i].rad)
      return 1;
  }
  return 0;
}

void Lmp_data::write_data(const char* outfilename){
  FILE *outfile;
  outfile = fopen(outfilename, "w");

  maxrad = -1;
  for (int i = 0; i < ntot; i++){
    if (particles[i].rad > maxrad) maxrad = particles[i].rad;
  }
  double dx, dy;
  dx = dy = 0;
  //if (!xperiodic) dx = maxrad;
  //if (!yperiodic) dy = maxrad;
  fprintf(outfile,"LAMMPS data file generated by place_particles\n");
  fprintf(outfile,"%d atoms\n\n",ntot);
  fprintf(outfile,"1 atom types\n\n");
  fprintf(outfile,"%g %g xlo xhi\n",xlo-dx,xhi+dx);
  fprintf(outfile,"%g %g ylo yhi\n",ylo-dy,yhi+dy);
  fprintf(outfile,"%g %g zlo zhi\n\n",zlo,zhi);

  fprintf(outfile,"Atoms\n\n");
  int id = 1;
  for (int i = 0; i < ntot; ++i){
        fprintf(outfile,"%d %d %g %g %g %g %g\n",
            id++, particles[i].type, particles[i].diam, particles[i].density,
            particles[i].x, particles[i].y, particles[i].z);
  }

  fclose(outfile);

}

