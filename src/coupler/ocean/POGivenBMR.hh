// Copyright (C) 2011 Constantine Khroulev
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

#ifndef _POGIVENBMR_H_
#define _POGIVENBMR_H_

#include "PGivenBMR.hh"
#include "POModifier.hh"

class POGivenBMR : public PGivenBMR<POModifier,PISMOceanModel>
{
public:
  POGivenBMR(IceGrid &g, const NCConfigVariable &conf)
    : PGivenBMR<POModifier,PISMOceanModel>(g, conf, NULL)
  {
    mass_flux_name  = "shelfbmassflux";
    option_prefix   = "-ocean_bmr";

    shelfbtemp.init_2d("shelfbtemp", g);
    shelfbtemp.set_string("pism_intent", "climate_state");
    shelfbtemp.set_string("long_name",
			  "absolute temperature at ice shelf base");
    shelfbtemp.set_units("Kelvin"); 

    shelfbmassflux.init_2d("shelfbmassflux", g);
    shelfbmassflux.set_string("pism_intent", "climate_state");
    shelfbmassflux.set_string("long_name",
    			      "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
    shelfbtemp.set_units("m s-1"); 
  }

  virtual ~POGivenBMR() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode sea_level_elevation(PetscReal &result) {
    result = sea_level;
    return 0;
  }

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);

  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,
  					  PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename); 

protected:
  string reference;
  bool adjust_bmr_set;
  PetscReal ref_openocean_shelfthk;
  IceModelVec2S *ice_thickness, ref_shelfbaseelev_array;
  NCSpatialVariable shelfbtemp, shelfbmassflux;
};


#endif /* _POGIVENBMR_H_ */
