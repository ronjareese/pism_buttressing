// Copyright (C) 2011, 2012 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _POGIVENTH_H_
#define _POGIVENTH_H_

#include "PGivenClimate.hh"
#include "POModifier.hh"

class POGivenTH : public PGivenClimate<POModifier,PISMOceanModel>
{
public:
  POGivenTH(IceGrid &g, const NCConfigVariable &conf)
    : PGivenClimate<POModifier,PISMOceanModel>(g, conf, NULL)
  {
    temp_name       = "thetao";
    mass_flux_name  = "salinity"; //NOTE: salinity_name instead of mass_flux_name
    option_prefix   = "-ocean_th";
  }

  virtual ~POGivenTH() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode sea_level_elevation(PetscReal &result) {
    result = sea_level;
    return 0;
  }

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);

  virtual PetscErrorCode calculate_boundlayer_temp_and_salt();

  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual PetscErrorCode shelf_base_temp_salinity_3eqn(vector<double> gat_array, PetscInt i, PetscInt j,
						       PetscReal gas, PetscReal sal_ocean,
						       PetscReal temp_insitu, PetscReal zice,
						       PetscReal &temp_base, PetscReal &sal_base);

  virtual PetscErrorCode compute_meltrate_3eqn(vector<double> gat_array, PetscReal gas, PetscReal rhow,
					       PetscReal rhoi, PetscReal temp_base, PetscReal sal_base,
					       PetscReal sal_ocean, PetscReal &meltrate);

  virtual PetscErrorCode adlprt(PetscReal salz, PetscReal temp_insitu, PetscReal pres, PetscReal &adlprt_out);
  virtual PetscErrorCode pttmpr(PetscReal salz, PetscReal temp_insitu, PetscReal pres,PetscReal rfpres, PetscReal &thetao);
  virtual PetscErrorCode potit(PetscReal salz,PetscReal thetao,PetscReal pres,PetscReal rfpres, PetscReal &temp_insitu_out);
  protected:
    IceModelVec2S *ice_thickness, *topg; // is not owned by this class
    IceModelVec2S temp_boundlayer, salinity_boundlayer;

    vector<double> gat_array;
    PetscReal gas;
    bool gamma_T_set;
    bool gamma_T_separate_set;
};

#endif /* _POGIVENTH_H_ */
