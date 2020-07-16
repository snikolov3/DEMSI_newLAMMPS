#include "demsi_remapping.h"
#include "demsi_communication.h"
#include "demsi_column.h"

//#define DEBUG_OBR

namespace DEMSI {

//------------------------------------------------------------------------------
// Remapping correction with OBR
//------------------------------------------------------------------------------

void Remapping::particles_in_use(const Mat Weights, DEMSI::ColumnVariable<double>* effectiveElementArea, std::vector<int> &useParticle) {

  // determine areas with initial weights
  double* effectiveAreaIce = new double[*(particles->nParticles)];
  for (int iParticleOld = 0 ; iParticleOld < *(particles->nParticles) ; iParticleOld++) {
    effectiveAreaIce[iParticleOld] = (*effectiveElementArea)(iParticleOld);
  } // iParticleInit

  double* effectiveAreaIceRemap = new double[tessellation->nParticlesAllInit];
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    effectiveAreaIceRemap[iParticleInit] = 0.0;
  } // iParticleInit

  remap_variable(1, effectiveAreaIceRemap, effectiveAreaIce, Weights);

  // find new elements in use
  useParticle.resize(tessellation->nParticlesAllInit);
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    if (effectiveAreaIceRemap[iParticleInit] > 0.0 and tessellation->typeInit[iParticleInit] == 0) {
      useParticle[iParticleInit] = 1;
    } else {
      useParticle[iParticleInit] = 0;
    }
  } // iParticleInit

} // Remapping::particles_in_use

// Correction area with OBR
void Remapping::correct_area_obr(DEMSI::ColumnVariable<double> *effectiveArea, std::vector<double> &extensiveAreaRemap, std::vector<double> &extensiveArea, std::vector<int> useParticle) {

  std::vector<double> areaMin;
  std::vector<double> areaMax;
  std::vector<double> areaVal;
  std::vector<double> parentVal;

  areaMin.resize(tessellation->nParticlesAllInit);
  areaMax.resize(tessellation->nParticlesAllInit);
  areaVal.resize(tessellation->nParticlesAllInit);
  parentVal.resize(tessellation->nParticlesAllInit);

  // check for particles on edge
  std::vector<int> partOnEdge;
  partOnEdge.resize(tessellation->nParticlesAllInit);
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    for (int iCellOnCell = 0 ; iCellOnCell < tessellation->nVerticesOnCellInit[iParticleInit] ; iCellOnCell++) {
      int iParticleInitNeighbour = tessellation->cellsOnCellLocal[iParticleInit][iCellOnCell];
     if (iParticleInitNeighbour != -1 and useParticle[iParticleInit] == 1 and useParticle[iParticleInitNeighbour] == 0) {
       partOnEdge[iParticleInit] = 1;
     }
    } // iCellOnCell
  } // iParticleInit

   int ipart = 0;
   double totalArea = 0;
   for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit; iParticleInit++) {
     parentVal[iParticleInit] = tessellation->initPolygons[iParticleInit].area();
     areaMax[iParticleInit] = 1.0;
     areaVal[iParticleInit] = extensiveAreaRemap[iParticleInit]/tessellation->initPolygons[iParticleInit].area();
     if (particles->type[iParticleInit] != 2){
      areaMin[iParticleInit] = std::min(areaVal[iParticleInit],1.0);
     }
     else {
       areaMin[iParticleInit] = 1.0;
     }
     totalArea += areaVal[iParticleInit]*parentVal[iParticleInit];
  } // iParticleInit


#ifdef DEBUG_OBR

  double Area = 0.0;
  for (int i = 0 ; i < extensiveArea.size() ; i++) {
    Area += extensiveArea[i];
  }

  std::cout << " Area_old = " << Area << " area_size = " << extensiveArea.size() << "\n";

  double Area0 = 0.0;
  for (int i = 0 ; i < areaVal.size() ; i++) {
    Area0 += areaMin[i]*tessellation->initPolygons[i].area();
  } // i
  std::cout << " Area_min = " << Area0 <<"\n";

  Area0 = 0.0;
  for (int i = 0 ; i < areaVal.size() ; i++) {
    Area0 += areaMax[i]*tessellation->initPolygons[i].area();
  } // i
  std::cout << " Area_max = " << Area0 <<"\n";

  Area0 = 0.0;
  for (int i = 0 ; i < areaVal.size() ; i++) {
      Area0 += areaVal[i]*tessellation->initPolygons[i].area();
  } // i
  std::cout << " Area_current = " << Area0 <<"\n";
#endif

  limit_obr_simple(areaVal, parentVal, areaMin, areaMax, totalArea);

  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    if (areaVal[iParticleInit] > 1e-10){
       extensiveAreaRemap[iParticleInit] = areaVal[iParticleInit]*parentVal[iParticleInit];
    }
    else {
       extensiveAreaRemap[iParticleInit] = 0.0;
    }
  }

#ifdef DEBUG_OBR
  Area = 0.0;
  for (int i = 0 ; i < extensiveAreaRemap.size() ; i++) {
      Area += extensiveAreaRemap[i];
  } // i

  std::cout << "Area_new = " << Area <<"\n";
#endif

} // Remapping::correct_area_obr

// Correction with OBR
void Remapping::correct_tracer_obr(DEMSI::ColumnVariable<double> *tracer, std::vector<double> &intensiveTracerRemap, std::vector<double> &extensiveTracerRemap, std::vector<double> &extensiveTracer, DEMSI::ColumnVariable<double> *parent, std::vector<double> &extensiveTracerParentRemap, std::vector<double> &extensiveTracerParent, std::vector<int> useParticle, std::vector<double> weights, std::vector<int> iParticlesOldWeights, std::vector<int> iParticlesInitWeights) {

  std::vector<double> tracerMin;
  std::vector<double> tracerMax;
  std::vector<double> tracerMinAgg;
  std::vector<double> tracerMaxAgg;
  std::vector<double> intTracerSum;

  // get local min/max of intensive tracer
  int sizePerParticle = tracer->size_per_particle();

  tracerMin.resize(tessellation->nParticlesAllInit*sizePerParticle);
  tracerMax.resize(tessellation->nParticlesAllInit*sizePerParticle);
  tracerMinAgg.resize(tessellation->nParticlesAllInit);
  tracerMaxAgg.resize(tessellation->nParticlesAllInit);
  intTracerSum.resize(tessellation->nParticlesAllInit);

  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    tracerMinAgg[iParticleInit] = 1e10;
    tracerMaxAgg[iParticleInit] = -1e10;
    for (int isizePerPart = 0 ; isizePerPart < sizePerParticle ; isizePerPart++) {
       int jind  = iParticleInit  * sizePerParticle + isizePerPart;
       tracerMin[jind] = 1e10;
       tracerMax[jind] = -1e10;
     } // sizePerPart
   } // iParticleInit

  int sizePerParticleParent = parent->size_per_particle();
  int sizeChildPerParent = tracer->size_per_particle() / parent->size_per_particle();

  // get aggregate old tracer value
  for (int iParticleOld = 0 ; iParticleOld < *(particles->nParticles) ; iParticleOld++) {
    intTracerSum[iParticleOld] = 0.0;
    for (int iParent = 0 ; iParent < sizePerParticleParent ; iParent++) {
      for (int iChild = 0 ; iChild < sizeChildPerParent ; iChild++) {
        int jOldc  = iParticleOld  * sizePerParticle + iParent * sizeChildPerParent + iChild;
        int jOldp  = iParticleOld  * sizePerParticleParent + iParent;
        if (extensiveTracerParent[jOldp] != 0) {
            intTracerSum[iParticleOld] += extensiveTracer[jOldc]/extensiveTracerParent[jOldp];
        }
      }
    }
  }

  // get min/max tracer value and aggregate values from old values - needs to be parallelized
  for (int iWeight = 0 ; iWeight < weights.size(); iWeight++) {

    int iParticleInit = iParticlesInitWeights[iWeight];
    int iParticleOld  = iParticlesOldWeights[iWeight];

    for (int iParent = 0 ; iParent < sizePerParticleParent ; iParent++) {
      for (int iChild = 0 ; iChild < sizeChildPerParent ; iChild++) {
        int jOldc  = iParticleOld  * sizePerParticle + iParent * sizeChildPerParent + iChild;
        int jOldp  = iParticleOld  * sizePerParticleParent + iParent;
        int jInit = iParticleInit * tracer->size_per_particle() + iParent * sizeChildPerParent + iChild;
  //     if (particles->type[iParticleOld] != 2){
        if (extensiveTracerParent[jOldp] != 0) {
           double tracerOld = extensiveTracer[jOldc]/extensiveTracerParent[jOldp];
           tracerMin[jInit] = std::min(tracerMin[jInit],tracerOld);
           tracerMax[jInit] = std::max(tracerMax[jInit],tracerOld);
        }
      } //iChild
    } //iParent

    tracerMinAgg[iParticleInit] = std::min(tracerMinAgg[iParticleInit],intTracerSum[iParticleOld]);
    tracerMaxAgg[iParticleInit] = std::max(tracerMaxAgg[iParticleInit],intTracerSum[iParticleOld]);

  } // iWeight

  // check for particles on edge of distribution
  std::vector<int> partOnEdge;
  partOnEdge.resize(tessellation->nParticlesAllInit);
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    partOnEdge[iParticleInit] = 0.0;
    for (int iCellOnCell = 0 ; iCellOnCell < tessellation->nVerticesOnCellInit[iParticleInit] ; iCellOnCell++) {
      int iParticleInitNeighbour = tessellation->cellsOnCellLocal[iParticleInit][iCellOnCell];
      if (iParticleInitNeighbour != -1 and useParticle[iParticleInit] == 1 and useParticle[iParticleInitNeighbour] == 0) {
        partOnEdge[iParticleInit] = 1;
      }
    } // iCellOnCell
  } // iParticleInit

#ifdef DEBUG_OBR
  if (tracer->name() == "iceAreaCategory") {
     double total_area_remap = 0;
     double total_area_init = 0;
     for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit; iParticleInit++) {
        total_area_remap = total_area_remap + extensiveTracerParentRemap[iParticleInit];
        if (extensiveTracerParentRemap[iParticleInit] > 0)
           total_area_init = total_area_init + tessellation->initPolygons[iParticleInit].area();
     } // iParticleInit
     std::cout << "\n" << "  total_area_remap = " << total_area_remap << "  total_area_init = " << total_area_init << "\n";
  }
#endif

 // Initialize arrays before setting final values
  std::vector<double> tracerVal;
  std::vector<double> parentVal;
  tracerVal.resize(tessellation->nParticlesAllInit*sizePerParticle);
  parentVal.resize(tessellation->nParticlesAllInit*sizePerParticle);
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit*sizePerParticle ; iParticleInit++) {
     parentVal[iParticleInit] = 0.0;
     tracerVal[iParticleInit] = 0.0;
     if (tracerMin[iParticleInit] > tracerMax[iParticleInit]) {
       tracerMin[iParticleInit] = 0.0;
       tracerMax[iParticleInit] = 0.0;
     }
  } // iParticleInit
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    if (tracerMinAgg[iParticleInit] > tracerMaxAgg[iParticleInit]) {
       tracerMinAgg[iParticleInit] = 0.0;
       tracerMaxAgg[iParticleInit] = 0.0;
    }
  }

 // ensure the aggregate  min/max of each area category is less than 1
 if (tracer->name() == "iceAreaCategory") {
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    double tracerMaxtest = 0.0;
    for (int isizePerPart = 0 ; isizePerPart < sizePerParticle ; isizePerPart++) {
      int jind  = iParticleInit  * sizePerParticle + isizePerPart;
      tracerMaxtest += tracerMax[jind];
    }
    if (tracerMaxtest > 1.0) {
       //std::cout << "iParticleInit = " << iParticleInit << " tracermax = " << tracerMaxtest << "\n";
       for (int isizePerPart = 0 ; isizePerPart < sizePerParticle ; isizePerPart++) {
         int jind  = iParticleInit  * sizePerParticle + isizePerPart;
           tracerMax[jind] = tracerMax[jind]*(1.0/tracerMaxtest);
           if (tracerMax[jind] < tracerMin[jind]){
              tracerMin[jind] = tracerMax[jind];
           }

             // Set min/max to remapped values in certain cases
             //tracerMax[jind] = intensiveTracerRemap[jind];
             //tracerMin[jind] = intensiveTracerRemap[jind];

//           std::cout << " tracerMax = " << tracerMax[jind] << "  tracerMin = " << tracerMin[jind] << "  tracerVal = " << intensiveTracerRemap[jind] << "\n";
       }
     } // isizePerPart
   } // iParticleInit
 } // if tracer = iceAreaCategory

/*
  // ensure the aggregate min/max of each category are satisfied for all tracers
  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    double tracerMintest = 0.0;
    double tracerMaxtest = 0.0;
    for (int isizePerPart = 0 ; isizePerPart < sizePerParticle ; isizePerPart++) {
      int jind  = iParticleInit  * sizePerParticle + isizePerPart;
      tracerMintest += tracerMin[jind];
      tracerMaxtest += tracerMax[jind];
    }
    if (tracerMintest < tracerMinAgg[iParticleInit]) {
      std::cout << "iParticleInit = " << iParticleInit << " tracermin = " << tracerMintest;
      std::cout << "   aggtracermin = " << tracerMinAgg[iParticleInit] << "\n";
      for (int isizePerPart = 0 ; isizePerPart < sizePerParticle ; isizePerPart++) {
        int jind  = iParticleInit  * sizePerParticle + isizePerPart;
        tracerMin[jind] = tracerMin[jind]*(tracerMinAgg[iParticleInit]/tracerMaxtest);
        if (tracerMax[jind] < tracerMin[jind]){
           tracerMax[jind] = tracerMin[jind];
        }
      }
    }
    if (tracerMaxtest > tracerMaxAgg[iParticleInit]) {
      std::cout << "iParticleInit = " << iParticleInit << " tracermax = " << tracerMaxtest;
      std::cout << "   aggtracermax = " << tracerMaxAgg[iParticleInit] << "\n";
      for (int isizePerPart = 0 ; isizePerPart < sizePerParticle ; isizePerPart++) {
        int jind  = iParticleInit  * sizePerParticle + isizePerPart;
        tracerMax[jind] = tracerMax[jind]*(tracerMaxAgg[iParticleInit]/tracerMaxtest);
        if (tracerMax[jind] < tracerMin[jind]){
          tracerMin[jind] = tracerMax[jind];
        }
      }
    }
  }
*/

#ifdef DEBUG_OBR
  std::cout << " tracer name = " << tracer->name()  << " parent name = " << parent->name() <<"\n";
  std::cout << " extensiveTracer size = " << extensiveTracer.size() << "  intensiveTracerRemap size = " << intensiveTracerRemap.size() <<"\n";
#endif

 for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
   for (int iParent = 0 ; iParent < sizePerParticleParent ; iParent++) {
     for (int iChild = 0 ; iChild < sizeChildPerParent ; iChild++) {
       int ic = iParticleInit * sizePerParticle + iParent * sizeChildPerParent + iChild;
       int ip  = iParticleInit  * sizePerParticleParent + iParent;

       tracerVal[ic] = intensiveTracerRemap[ic];

       // for tracers/flux whose parent value is effectiveElementArea, set to initial polygon area
       double ratio = 1;
       if (parent->name() == "effectiveElementArea") {
         parentVal[ic] = tessellation->initPolygons[iParticleInit].area();
         ratio = extensiveTracerParentRemap[iParticleInit]/tessellation->initPolygons[iParticleInit].area();
         if (partOnEdge[iParticleInit] == 1) {
            ratio = std::min(ratio,1.0);
            tracerVal[ic] = tracerVal[ic]*ratio;
         }
       }
       else {
         parentVal[ic] = extensiveTracerParentRemap[ip];
       }

       // if particle is on edge of distribution set min/max to 0 depending on sign of value
       if (partOnEdge[iParticleInit] == 1) {
         tracerMin[ic]  = std::min(tracerMin[ic],0.0);
         tracerMax[ic]  = std::max(tracerMax[ic],0.0);
       }
     } // iChild
   } // iParent
 } // iParticleInit

  double Mass = 0.0;
  for (int i = 0 ; i < extensiveTracer.size() ; i++) {
      Mass += extensiveTracer[i];
  }

#ifdef DEBUG_OBR
  std::cout << " Total_Tracer_Mass_old = " << Mass << " tracer_size = " << extensiveTracer.size() << "\n";

  double Mass0 = 0.0;
  for (int i = 0 ; i < tracerVal.size() ; i++) {
      Mass0 += tracerMin[i]*parentVal[i];
  } // i
  std::cout << " Total_Tracer_Mass_min = " << Mass0 <<"\n";

  Mass0 = 0.0;
  for (int i = 0 ; i < tracerVal.size() ; i++) {
      Mass0 += tracerMax[i]*parentVal[i];
  } // i
  std::cout << " Total_Tracer_Mass_max = " << Mass0 <<"\n";

  Mass0 = 0.0;
  for (int i = 0 ; i < tracerVal.size() ; i++) {
      Mass0 += tracerVal[i]*parentVal[i];
  } // i
  std::cout << " Total_Tracer_Mass_current = " << Mass0 <<"\n";

  Mass0 = 0.0;
  for (int i = 0 ; i < extensiveTracerRemap.size() ; i++) {
      Mass0 += extensiveTracerRemap[i];
  } // i
  std::cout << " Total_Tracer_Mass_remap = " << Mass0 <<"\n";

#endif

  //limit_obr(tracerVal, parentVal, tracerMin, tracerMax, Mass);
  limit_obr_simple(tracerVal, parentVal, tracerMin, tracerMax, Mass);

  for (int iParticleInit = 0 ; iParticleInit < tessellation->nParticlesAllInit ; iParticleInit++) {
    for (int i = 0 ; i < sizePerParticle ; i++) {
      int j = iParticleInit * sizePerParticle + i;
         extensiveTracerRemap[j] = tracerVal[j]*parentVal[j];
    } // i
  } // iParticleInit

#ifdef DEBUG_OBR
  Mass = 0.0;
  for (int i = 0 ; i < extensiveTracerRemap.size() ; i++) {
      Mass += extensiveTracerRemap[i];
  } // i

  std::cout << " Total_Tracer_Mass_new = " << Mass <<"\n";
#endif

} // Remapping::correct_tracer_obr

//  simplified version of optimization based remap algorithm
void Remapping::limit_obr_simple(std::vector<double> &tracerVal, std::vector<double> &parentVal, std::vector<double> &tracerMin, std::vector<double> &tracerMax, double Mass){

  int maxclip = 50;
  int nclip = 0;
  double eta = 1e-11;
  double hfd = 1e-06;
  double lambda = 0;
  double rp,rc,rd,lambda_p,lambda_c,alpha;

  if (std::abs(Mass) < 1e-13) return;

  // create a temporary array to store new column variable
  std::vector<double> tmp;
  tmp.resize(tracerVal.size());

  // get residual
  rp = 0;
  solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, rp);
  if (std::abs(rp) < eta) {
     std::cout << "Zero residual in OBR: " << rp << "\n";
     return;
  }

  //increment counter
  nclip = nclip + 1;

  // get residual
  rc = 0;
  lambda = lambda + hfd;
  solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, rc);

  // increment counter
  nclip = nclip + 1;

  // finite difference slope
  if (std::abs(rc-rp) < 1e-15) {
    if (std::abs(rp) > eta) {
       std::cout << "Nonzero residual in OBR: " << rc << " , " << rp << "  nclip = " << nclip << "\n";
    }
    for (int i = 0 ; i < tracerVal.size() ; i++) {
      tracerVal[i] = tmp[i];
    }
    return;
  }
  alpha = hfd/(rc-rp);

  // finite difference approximation
  lambda_p = 0;
  lambda_c = lambda_p - alpha*rp;

  while(std::abs(rc) > eta && nclip < maxclip) {

    // get residual
    rc = 0;
    lambda = lambda_c;
    solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, rc);

    nclip = nclip + 1;

    if (std::abs(rp-rc) < 1e-15) {
      if (std::abs(rc) > eta) {
        std::cout << "Nonzero residual in OBR: " << rc << " , " << rp << "  nclip = " << nclip << "\n";
      }
      for (int i = 0 ; i < tracerVal.size() ; i++) {
        tracerVal[i] = tmp[i];
      }
      return;
    }
    alpha = (lambda_p - lambda_c)/(rp - rc);
    rp = rc;

    // take step
    lambda_p = lambda_c;
    lambda_c = lambda_c - alpha*rc;
    std::cout << "  nclip  = " << nclip << " rp = " << rp << "  alpha = " << alpha << " lambda_p = " << lambda_p << "\n";
  }

  for (int i = 0 ; i < tracerVal.size() ; i++) {
    tracerVal[i] = tmp[i];
  }
  std::cout << "  obr finish " <<"\n";

} // Remapping::limit_obr_simple

// optimization based remap, based on D. Ridzal implementation of
// Dai/Fletcher algorithm for singly linearly constrained quadratic
// programs with simple bound constraints
void Remapping::limit_obr(std::vector<double> &tracerVal, std::vector<double> &parentVal, std::vector<double> &tracerMin, std::vector<double> &tracerMax, double Mass){

  int maxclip = 30;
  int nclip = 0;
  double eta = 1e-11;
  double eta2 = 1e-11;
  double lambda = 0;
  double lambda_l, lambda_u, r_l, r_u, lambda_new, s;
  double Dlambda = 2;

   if (abs(Mass) < 1e-13) return;

   // create a temporary array to store new column variable
   std::vector<double> tmp;
   tmp.resize(tracerVal.size());

   // bracketing phase
   // get residual
   double r = 0;
   solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, r);
   std::cout << "Initial residual in OBR: " << r << "\n";
   if (std::abs(r) < eta) {
      std::cout << "Zero residual in OBR: " << r << "\n";
      return;
   }

   if (r < -eta){
     lambda_l = lambda;
     r_l = r;
     lambda = lambda + Dlambda;
     solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, r);
     nclip = nclip + 1;
     int iclip = 0;
     while (r < -eta2 && iclip < 10) {
       iclip = iclip+1;
       lambda_l = lambda;
       r_l = r;
       s = std::max(r_l/r-1.0,0.1);
       Dlambda = Dlambda + Dlambda/s;
       lambda = lambda + Dlambda;
       solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, r);
       nclip = nclip + 1;
     }
     lambda_u = lambda;
     r_u = r;
     std::cout << " Bracketing:  r_u = " << r_u << "  lambda_u = " << lambda_u << "\n";
     std::cout << " Bracketing:  r_l = " << r_l << "  lambda_l = " << lambda_l << "\n";
   }
   else if (r > eta) {
     lambda_u = lambda;
     r_u = r;
     lambda = lambda - Dlambda;
     solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, r);
     nclip = nclip + 1;
     int iclip = 0;
     while (r > eta2 && iclip < 10) {
       iclip = iclip+1;
       lambda_u = lambda;
       r_u = r;
       s = std::max(r_u/r-1.0,0.1);
       Dlambda = Dlambda + Dlambda/s;
       lambda = lambda - Dlambda;
       solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, r);
       nclip = nclip + 1;
     }
     lambda_l = lambda;
     r_l = r;
     std::cout << " Bracketing:  r_l = " << r_l << "  lambda_l = " << lambda_l << "\n";
     std::cout << " Bracketing:  r_u = " << r_u << "  lambda_u = " << lambda_u << "\n";
   }
   else {
     for (int i = 0 ; i < tracerVal.size() ; i++) {
       tracerVal[i] = tmp[i];
     }
     return;
   }

   //secant phase
   s = 1.0 - r_l/r_u;
   Dlambda = Dlambda/s;
   lambda = lambda_u - Dlambda;

   // get residual
   solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, r);

   // increment counter
   nclip = nclip + 1;

   while(std::abs(r) > eta && nclip < maxclip) {

     if (r > eta2) {
       if (s <= 2.0) {
          lambda_u = lambda;
          r_u = r;
          s = 1.0-r_l/r_u;
          Dlambda = (lambda_u-lambda_l)/s;
          lambda = lambda_u - Dlambda;
       }
       else {
          s = std::max(r_u/r - 1.0,0.1);
          Dlambda = (lambda_u-lambda)/s;
          lambda_new = std::max(lambda-Dlambda, 0.75*lambda_l+0.25*lambda);
          lambda_u = lambda;
          r_u = r;
          lambda = lambda_new;
          s = (lambda_u-lambda_l)/(lambda_u-lambda);
       }
     }
     else if (r < -eta2) {
        if (s >= 2.0) {
          lambda_l = lambda;
          r_l = r;
          s = 1.0 - r_l/r_u;
          Dlambda = (lambda_u-lambda_l)/s;
          lambda = lambda_u - Dlambda;
        }
        else {
          s = std::max(r_l/r-1.0, 0.1);
          Dlambda = (lambda - lambda_l)/s;
          lambda_new = std::min(lambda+Dlambda, 0.75*lambda_u+0.25*lambda);
          lambda_l = lambda;
          r_l = r;
          lambda = lambda_new;
          s = (lambda_u-lambda_l)/(lambda_u-lambda);
        }
      }
      else {
        for (int i = 0 ; i < tracerVal.size() ; i++) {
          tracerVal[i] = tmp[i];
        }
        return;
      }
      solve_qp(tmp, tracerVal, parentVal, tracerMin, tracerMax, Mass, lambda, r);

      nclip = nclip + 1;

      std::cout << "  nclip  = " << nclip << " r = " << r << "  lambda = " << lambda << "\n";
  }

  for (int i = 0 ; i < tracerVal.size() ; i++) {
     tracerVal[i] = tmp[i];
  }
     std::cout << "  obr finish " <<"\n";

} // Remapping::limit_obr

// solve quadratic program
void Remapping::solve_qp(std::vector<double> &tracerArray, std::vector<double> &tracerVal, std::vector<double> &parentVal, std::vector<double> &tracerMin,
                         std::vector<double> &tracerMax, double Mass, double lambda, double &residual){


   for (int i = 0 ; i < tracerArray.size(); i++) {
      tracerArray[i] = tracerVal[i] + lambda*parentVal[i]/Mass;
   }

  // limit - take median
   for (int i = 0 ; i < tracerArray.size(); i++) {
     if (tracerMin[i] > tracerMax[i]) {
        tracerArray[i] = 0;
     }
     else {
        tracerArray[i] = std::min(std::max(tracerArray[i], tracerMin[i]), tracerMax[i]);
     }
   }

  // compute residual  (scaled by mass) - need to parallelize
   residual = 0;
   for (int i = 0 ; i < tracerArray.size(); i++) {
      residual += tracerArray[i]*parentVal[i]/Mass;
   }
   residual = residual - 1;

} // Remapping::solve_qp

} // namespace DEMSI
