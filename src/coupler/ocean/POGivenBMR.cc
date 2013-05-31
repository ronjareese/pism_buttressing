// Copyright (C) 2011, 2012 Constantine Khroulev
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

#include "POGivenBMR.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"
#include "pism_options.hh" 

PetscErrorCode POGivenBMR::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool regrid = true;
  int start = -1;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the ocean model 'BMR' reading sub-shelf mass flux from a file...\n"); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters(""); CHKERRQ(ierr);

  // ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr);

  // ierr = temp.set_attrs("climate_forcing",
  //                       "absolute temperature at ice shelf base",
  //                       "Kelvin", ""); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing",
                       "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                       "m s-1", ""); CHKERRQ(ierr);

  // ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);


  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(grid.com, 1, "ERROR: ice thickness is not available"); }

  // read time-independent data right away:
  // if (temp.get_n_records() == 1 && mass_flux.get_n_records() == 1) {
  //   ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  // }

  // read time-independent data right away:
  if (mass_flux.get_n_records() == 1) {
    ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  }


  // RTOPO-1
  ierr = PISMOptionsIsSet("-adjust_bmr", adjust_bmr_set); CHKERRQ(ierr);
  if (adjust_bmr_set) { 
    ierr = verbPrintf(2, grid.com,
                      "* Sub-shelf mass flux will be adjusted according to reference ice shelf base elevation"); CHKERRQ(ierr);   

    ierr = find_pism_input(filename, regrid, start); CHKERRQ(ierr);

    ierr = RTOPO1_shelfbaseelev.create(grid, "draft", false); CHKERRQ(ierr);
    ierr = RTOPO1_shelfbaseelev.set_attrs("climate_state",
					  "rtopo1 ice geometry",
					  "m",
					  ""); CHKERRQ(ierr); // no CF standard_name ??
    // ierr = RTOPO1_shelfbaseelev.set_glaciological_units("m");
    // RTOPO1_shelfbaseelev.write_in_glaciological_units = true;
    // RTOPO1_shelfbaseelev.time_independent = true;

    // ierr = find_pism_input(rtopo1_filename, regrid, start); CHKERRQ(ierr);
    
    // read RTOPO-1 ice geometry from file
    ierr = verbPrintf(2, grid.com, 
		      "\n  - Reading RTopo-1 data from '%s' ... \n",
		      filename.c_str()); CHKERRQ(ierr); 
    if (regrid) {
      ierr = RTOPO1_shelfbaseelev.regrid(filename.c_str(), true); CHKERRQ(ierr); // fails if not found!
    } else {
      ierr = RTOPO1_shelfbaseelev.read(filename.c_str(), start); CHKERRQ(ierr); // fails if not found!
    }
    string rtopo1_history = "read from " + filename + "\n";

    ierr = RTOPO1_shelfbaseelev.set_attr("history", rtopo1_history); CHKERRQ(ierr);
  }

  // delete lic;

  return 0;
}

PetscErrorCode POGivenBMR::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = mass_flux.at_time(t); CHKERRQ(ierr);
  // ierr = temp.at_time(t); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POGivenBMR::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const PetscScalar T0 = config.get("water_melting_point_temperature"), // K
    beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity"), // K m-1
    ice_rho = config.get("ice_density"),
    sea_water_rho = config.get("sea_water_density");

  // ierr = ice_thickness->begin_access();   CHKERRQ(ierr);
  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // const PetscScalar shelfbaseelev = - ( ice_rho / sea_water_rho ) * (*ice_thickness)(i,j); // FIXME issue #15
      const PetscScalar shelfbaseelev = - ( ice_rho / sea_water_rho ) * H[i][j]; // FIXME issue #15
      // temp is set to melting point at depth
      result(i,j) = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
    }
  }
  ierr = ice_thickness->end_access();   CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);  

  // PetscErrorCode ierr = temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode POGivenBMR::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const PetscScalar beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity"), // K m-1
    ice_rho = config.get("ice_density"),
    sea_water_rho = config.get("sea_water_density");

  PetscScalar **H;
  PetscReal dT_melt;
  PetscReal C_BMR = 10.0; // m/(a*K)
  PetscReal dH;

  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = mass_flux.begin_access(); CHKERRQ(ierr);
  if (adjust_bmr_set) ierr = RTOPO1_shelfbaseelev.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar PISM_shelfbaseelev = - ( ice_rho / sea_water_rho ) * H[i][j]; // FIXME issue #15
      // temp is set to melting point at depth
      // T_PISMmelt = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
      if (adjust_bmr_set) {
	dH = ( RTOPO1_shelfbaseelev(i,j) - PISM_shelfbaseelev );
	dT_melt = beta_CC_grad * ( RTOPO1_shelfbaseelev(i,j) - PISM_shelfbaseelev );
	// if (i==70 && j==100) {
	//   ierr = verbPrintf(2, grid.com, "\n dH = %f", dH); CHKERRQ(ierr);
	//   ierr = verbPrintf(2, grid.com, "\n PISM_shelfbaseelev = %f", PISM_shelfbaseelev); CHKERRQ(ierr);
	//   ierr = verbPrintf(2, grid.com, "\n RTOPO1_shelfbaseelev = %f", RTOPO1_shelfbaseelev(i,j)); CHKERRQ(ierr);
	//   ierr = verbPrintf(2, grid.com, "\n dT_melt = %f", dT_melt); CHKERRQ(ierr);
	// }
	// difference in shelf base elevation defines difference in pressure melting point
	result(i,j) = mass_flux(i,j) + ( C_BMR * dT_melt ) / secpera;
      } else {
	result(i,j) = mass_flux(i,j);
      }
      // dT_melt = 1.0;
      // dT_melt = T_RTopo1melt(i,j)-T_PISMmelt;
      // mass_flux(i,j) = mass_flux(i,j) + C_BMR * dT_melt;

      // result(i,j) = mass_flux(i,j) + 1.0/secpera; 
      // result(i,j) = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = mass_flux.end_access(); CHKERRQ(ierr);
  if (adjust_bmr_set) ierr = RTOPO1_shelfbaseelev.begin_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}

// void POGivenBMR::add_vars_to_output(string, map<string,IceModelVec2S> &result) {
void POGivenBMR::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {
  if (keyword == "medium" || keyword == "big") {
  result["shelfbtemp"] = shelfbtemp;
  // result["shelfbmassflux"] = shelfbmassflux;
  // result["shelfbmassflux"] = mass_flux;
  }
}

PetscErrorCode POGivenBMR::define_variables(set<string> vars, const PIO &nc,
                                            PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "shelfbtemp")) {
    ierr = shelfbtemp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  // if (set_contains(vars, "shelfbmassflux")) {
  //   // ierr = mass_flux.define(nc, nctype, true); CHKERRQ(ierr);
  //   ierr = shelfbmassflux.define(nc, nctype, true); CHKERRQ(ierr);
  // }

  return 0;
}

PetscErrorCode POGivenBMR::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbtemp, 0); CHKERRQ(ierr);
    ierr = shelf_base_temperature(tmp); CHKERRQ(ierr);
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  // if (set_contains(vars, "shelfbmassflux")) {
  //   if (!tmp.was_created()) {
  //     ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
  //   }

  //   ierr = tmp.set_metadata(mass_flux, 0); CHKERRQ(ierr);
  //   ierr = tmp.set_metadata(shelfbmassflux, 0); CHKERRQ(ierr);
  //   tmp.write_in_glaciological_units = true;
  //   ierr = shelf_base_mass_flux(tmp); CHKERRQ(ierr);
  //   ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  // }

  return 0;
}
