#include "demsi_column.h"

#include "demsi_logging.h"
#include "demsi_particles.h"

#include <iomanip>

extern "C" {void ridging_set_variables(int*, double*, double*, double*, double*, double*, double*);}
extern "C" {void ridging_ridge(int*, double*, double*, double*, double*, double*, double*, double*);}

namespace DEMSI {

// set ridging variables needed for dynamics
void Column::set_ridging_variables(void) {

  if (useColumn and useColumnRidging) {

    int nCategories = columnDimensions->size("nCategories");

    // particles variables on device
    auto ridgingIceThickness_particles        = particles->ridgingIceThickness;
    auto ridgingIceThicknessWeight_particles  = particles->ridgingIceThicknessWeight;
    auto netToGrossClosingRatio_particles     = particles->netToGrossClosingRatio;
    auto changeEffectiveElementArea_particles = particles->changeEffectiveElementArea;

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      double ridgingIceThickness;
      double ridgingIceThicknessWeight;
      double netToGrossClosingRatio;

      ridging_set_variables(&nCategories,
			    openWaterArea->get(iParticle),
			    iceAreaCategory->get(iParticle),
			    iceVolumeCategory->get(iParticle),
			    &ridgingIceThickness,
			    &ridgingIceThicknessWeight,
			    &netToGrossClosingRatio);

      ridgingIceThickness_particles       (iParticle) = ridgingIceThickness;
      ridgingIceThicknessWeight_particles (iParticle) = ridgingIceThicknessWeight;
      netToGrossClosingRatio_particles    (iParticle) = netToGrossClosingRatio;
      changeEffectiveElementArea_particles(iParticle) = 0.0;

    } // iParticle

  } // useColumnRidging

} // Column::set_ridging_variables

// perform column ridging
void Column::ridge(void) {

  if (useColumn and useColumnRidging) {

    int nCategories = columnDimensions->size("nCategories");

    auto netToGrossClosingRatio_particles     = particles->netToGrossClosingRatio;
    auto changeEffectiveElementArea_particles = particles->changeEffectiveElementArea;

    for (int iParticle = 0 ; iParticle < *(particles->nParticles) ; iParticle++) {

      double netToGrossClosingRatio     = netToGrossClosingRatio_particles(iParticle);
      double changeEffectiveElementArea = changeEffectiveElementArea_particles(iParticle);

      ridging_ridge(&nCategories,
		    openWaterArea->get(iParticle),
		    iceAreaCategory->get(iParticle),
		    iceVolumeCategory->get(iParticle),
		    categoryThicknessLimits,
		    &netToGrossClosingRatio,
		    effectiveElementArea->get(iParticle),
		    &changeEffectiveElementArea);

    } // iParticle

    // aggregate tracers
    aggregate();

  } // useColumnRidging

} // Column::ridge

} // namespace DEMSI
