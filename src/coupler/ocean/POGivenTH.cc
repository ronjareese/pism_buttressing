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

#include "POGivenTH.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"

/* TO DO:
 *
 * Translate all comments into English.
 *
 * Add detailed references ("this implements equations (3)-(7) in Smith and
 * Jones, 2001" and similar).
 *
 * Give meaningful names to class methods (potit? pttmpr? adlprt?). Long names
 * are OK.
 *
 * All computationally expensive code should be called from the update() method.
 *
 * Make sure that the code compiles without warnings.
 *
 * Remove unused code. (It can always be recovered from an earlier version of
 * the code.)
 */

PetscErrorCode POGivenTH::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the 3eqn melting parameterization ocean model\n"
                    "  reading ocean temperature and salinity from a file...\n"); CHKERRQ(ierr);

  gat_array.resize(1);
  gat_array[0] = config.get("gamma_T") ;  // give gat_array[0] a value (default or prescribed)
  gas = config.get("gamma_S");

  // if option -gamma_T is set, print value of gamma_T
  ierr = PISMOptionsIsSet("-gamma_T", gamma_T_set); CHKERRQ(ierr);
  if (gamma_T_set) {
    ierr = verbPrintf(2, grid.com,
                      "* Turbulent heat exchange coefficient for 3eqn melting parameterization\n"
                      "  is set for whole domain to gamma_T=%f \n", gat_array[0]); CHKERRQ(ierr);
  }

  // if option -gamma_T_separate is set, print list of gamma_T values
  ierr = PISMOptionsIsSet("-gamma_T_separate", gamma_T_separate_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-gamma_T", gamma_T_set); CHKERRQ(ierr);

  if (gamma_T_separate_set) {
    gat_array.resize(3);
    // default values:
    gat_array[0] = config.get("gamma_T");
    gat_array[1] = config.get("gamma_T");
    gat_array[2] = config.get("gamma_T");

    ierr = PISMOptionsRealArray("-gamma_T_separate", "gamma_T_PIG_North, gamma_T_PIG_South, gamma_T_TG",
                              gat_array, gamma_T_separate_set); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com,
                      "* Turbulent heat exchange coefficients for 3eqn melting parameterization\n"
                      "  are set separately to\n"
                      "  gamma_T_PIG_North=%f, gamma_T_PIG_South=%f, gamma_T_TG=%f \n", gat_array[0], gat_array[1], gat_array[2]); CHKERRQ(ierr);

    if (gat_array.size() != 3) {
      PetscPrintf(grid.com,
                "PISM ERROR: option -gamma_T_separate requires a comma-separated list with 3 numbers; got %d\n",
                gat_array.size());
      PISMEnd();
    }
  }

  if (gamma_T_set && gamma_T_separate_set) {
    PetscPrintf(grid.com,
                "PISM ERROR: Options -gamma_T AND -gamma_T_separate are used at the same time,\n"
                "this could cause trouble... aborting\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  // print value of gamma_S
  ierr = verbPrintf(2, grid.com, "gamma_S=%e \n", gas); CHKERRQ(ierr);

  ierr = process_options(); CHKERRQ(ierr);

  ierr = set_vec_parameters("", ""); CHKERRQ(ierr);

  ierr = temp.create(grid, temp_name, false); CHKERRQ(ierr);
  ierr = mass_flux.create(grid, mass_flux_name, false); CHKERRQ(ierr); //NOTE: salinity instead of mass_flux

  ierr = salinity_boundlayer.create(grid, "salt_bl", false); CHKERRQ(ierr);
  ierr = temp_boundlayer.create(grid, "temp_bl", false); CHKERRQ(ierr);

  ierr = temp.set_attrs("climate_forcing", "absolute temperature at ice shelf base",
                        "Kelvin", ""); CHKERRQ(ierr); // temp_boundlayer is not a standard name
  ierr = mass_flux.set_attrs("climate_forcing", "ocean salinity at ice shelf",
                             "g/kg",
                             ""); CHKERRQ(ierr); // salinity_boundlayer is not a standard name

  ierr = temp.init(filename); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) {SETERRQ(grid.com, 1, "ERROR: ice thickness is not available");}

  topg = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (!topg) {SETERRQ(grid.com, 1, "ERROR: bedrock altitude is not available");}

  // read time-independent data right away:
  if (temp.get_n_records() == 1 && mass_flux.get_n_records() == 1) {
    ierr = update(grid.time->current(), 0); CHKERRQ(ierr); // dt is irrelevant
  }


  // modified by Ricarda
  ierr = PISMOptionsIsSet("-ocean_th_deltaT", ocean_th_deltaT_set); CHKERRQ(ierr);
  
  if (ocean_th_deltaT_set) {
      bool delta_T_set;
      //bool delta_T_factor_set;
      string delta_T_file;
      //PetscScalar delta_T_factor;
  
      ierr = PISMOptionsString("-ocean_th_deltaT",
                               "Specifies the ocean temperature offsets file to use with -ocean_th_deltaT",
                               delta_T_file, delta_T_set); CHKERRQ(ierr);
                               
      //ierr = PISMOptionsReal(PETSC_NULL,"-ocean_th_factor",delta_T_factor, delta_T_factor_set); CHKERRQ(ierr);
                                                        
  
      ierr = verbPrintf(2, grid.com, 
                        "  reading delta_T data from forcing file %s for -ocean_th_deltaT actions ...\n",
                        delta_T_file.c_str());  CHKERRQ(ierr);
  
      delta_T = new Timeseries(grid.com, grid.rank, "delta_T",
                               grid.config.get_string("time_dimension_name"));
      ierr = delta_T->set_units("Kelvin", ""); CHKERRQ(ierr);
      ierr = delta_T->set_dimension_units(grid.time->units(), ""); CHKERRQ(ierr);
      ierr = delta_T->set_attr("long_name", "ocean temperature offsets");
      CHKERRQ(ierr);
      ierr = delta_T->read(delta_T_file, grid.time->use_reference_date()); CHKERRQ(ierr);
  }
  



  return 0;
}

PetscErrorCode POGivenTH::update(PetscReal my_t, PetscReal my_dt) {

  PetscErrorCode ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  ierr = mass_flux.at_time(t+0.5*dt); CHKERRQ(ierr);
  ierr = temp.at_time(t+0.5*dt); CHKERRQ(ierr);

  ierr = calculate_boundlayer_temp_and_salt(); CHKERRQ(ierr);

  //  ierr = mass_flux.average(t, dt); CHKERRQ(ierr);
  //  ierr = temp.average(t, dt); CHKERRQ(ierr);

  //  if (enable_time_averaging) {
  //    ierr = mass_flux.average(t, dt); CHKERRQ(ierr);
  //    ierr = temp.average(t, dt); CHKERRQ(ierr);
  //  } else {
  //    ierr = mass_flux.get_record_years(t); CHKERRQ(ierr);
  //    ierr = temp.get_record_years(t); CHKERRQ(ierr);
  //  }

  return 0;
}

//NOTE Ported from Matthias
//NOTE is this needed?
//PetscErrorCode POGivenTH::update(PetscReal t_years, PetscReal dt_years) {
//  PetscErrorCode ierr = update_internal(t_years, dt_years); CHKERRQ(ierr);
//
//  if (enable_time_averaging) {
//    ierr = mass_flux.average(t, dt); CHKERRQ(ierr);
//    ierr = temp.average(t, dt); CHKERRQ(ierr);
//  } else {
//    ierr = mass_flux.get_record_years(t); CHKERRQ(ierr);
//    ierr = temp.get_record_years(t); CHKERRQ(ierr);
//  }

//  return 0;
//}

//NOTE Ported from Matthias
PetscErrorCode POGivenTH::shelf_base_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr = temp_boundlayer.copy_to(result); CHKERRQ(ierr);
  return 0;
}

//NOTE Ported from Matthias
PetscErrorCode POGivenTH::calculate_boundlayer_temp_and_salt() {

  PetscErrorCode ierr;

  const PetscScalar rhoi = config.get("ice_density");
  const PetscScalar rhow = config.get("sea_water_density");
  const PetscScalar reference_pressure = 1.01325; // pressure of atmosphere in bar

  PetscReal temp_insitu, temp_base, sal_base;

  ierr = ice_thickness->begin_access();   CHKERRQ(ierr);
  ierr = topg->begin_access();   CHKERRQ(ierr);
  ierr = temp.begin_access(); CHKERRQ(ierr);
  ierr = mass_flux.begin_access(); CHKERRQ(ierr); // NOTE salinity instead of mass_flux
  ierr = salinity_boundlayer.begin_access(); CHKERRQ(ierr);
  ierr = temp_boundlayer.begin_access(); CHKERRQ(ierr);
  
  bool delta_T_factor_set;
  PetscReal delta_T_factor=1.0;
  if (ocean_th_deltaT_set) {
    ierr = PISMOptionsReal("-ocean_th_factor","-ocean_th_factor",delta_T_factor, delta_T_factor_set); CHKERRQ(ierr);
  }
  

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar shelfbaseelev = -(rhoi/rhow) * (*ice_thickness)(i,j); // FIXME issue #15
      
      // modified by Ricarda for the thetao output
      if ((delta_T != NULL) && ocean_th_deltaT_set) {
        temp(i,j) = temp(i,j) + (*delta_T)(t + 0.5 * dt) / delta_T_factor;
      }    
      
      PetscReal thetao    = temp(i,j) - 273.15; // to degC
      PetscReal sal_ocean = mass_flux(i,j); //NOTE salinity instead of mass_flux
      PetscReal press = rhow * -1. * shelfbaseelev/1000. + reference_pressure; //NOTE Unit??
      potit(sal_ocean, thetao, press, reference_pressure, temp_insitu);


      // ensure that melting at grounded cells do not yield too large pressure
      PetscReal zice = -1*PetscMin(PetscAbs(shelfbaseelev),PetscAbs((*topg)(i,j)));
      //ierr = verbPrintf(2, grid.com, "zice=%e \n", zice); CHKERRQ(ierr);

      //shelf_base_temp_salinity_3eqn(gat, sal_ocean, temp_insitu, zice,
      //                                temp_base, sal_base);

      shelf_base_temp_salinity_3eqn(gat_array, i, j, gas, sal_ocean, temp_insitu, zice, temp_base,
                                    sal_base);

      temp_boundlayer(i,j)     = temp_base + 273.15; // to Kelvin
      salinity_boundlayer(i,j) = sal_base;

    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr); // NOTE is this way ok? (Yes. -- CK)
  ierr = topg->end_access(); CHKERRQ(ierr);
  ierr = temp.end_access(); CHKERRQ(ierr);
  ierr = mass_flux.end_access(); CHKERRQ(ierr);
  ierr = salinity_boundlayer.end_access(); CHKERRQ(ierr);
  ierr = temp_boundlayer.end_access(); CHKERRQ(ierr);

  return 0;
}

//NOTE This is replaced (see below), copy_to ???
//PetscErrorCode POGivenTH::shelf_base_salinity(IceModelVec2S &result) {
//  PetscErrorCode ierr = salinity.copy_to(result); CHKERRQ(ierr);
//  return 0;
//}

//NOTE Ported from Matthias
PetscErrorCode POGivenTH::shelf_base_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  //   PetscReal L = config.get("water_latent_heat_fusion"),
  const PetscScalar rhow = config.get("sea_water_density");
  const PetscScalar rhoi = config.get("ice_density");
  //const PetscScalar gat = config.get("gamma_T");   // set gat from option "-gamma_T"
  PetscReal meltrate_3eqn, sal_ocean, sal_base, temp_base;

  PetscReal reference_pressure = 1.01325; // pressure of atmosphere in bar

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = mass_flux.begin_access(); CHKERRQ(ierr);//NOTE salinity instead of mass_flux
  ierr = salinity_boundlayer.begin_access(); CHKERRQ(ierr);
  ierr = temp_boundlayer.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      sal_ocean = mass_flux(i,j); //NOTE salinity instead of mass_flux
      sal_base  = salinity_boundlayer(i,j);
      temp_base = temp_boundlayer(i,j) - 273.15; // to degC

      compute_meltrate_3eqn(gat_array, gas, rhow, rhoi, temp_base, sal_base, sal_ocean, meltrate_3eqn);
      //mpute_meltrate_3eqn(rhow, rhoi, temp_base, sal_base, sal_ocean, meltrate_3eqn);
      result(i,j) = -1.0*meltrate_3eqn;

    }
  }
  //  ierr = mass_flux.copy_to(result); CHKERRQ(ierr);

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = mass_flux.end_access(); CHKERRQ(ierr); //NOTE salinity instead of mass_flux
  ierr = salinity_boundlayer.end_access(); CHKERRQ(ierr);
  ierr = temp_boundlayer.end_access(); CHKERRQ(ierr);

  return 0;
}

// Ported from Matthias
PetscErrorCode POGivenTH::shelf_base_temp_salinity_3eqn(vector<double> gat_array, PetscInt i,
							PetscInt j, PetscReal gas, PetscReal sal_ocean,
                                                        PetscReal temp_insitu, PetscReal zice,
							PetscReal &temp_base, PetscReal &sal_base){
//PetscErrorCode POGivenTH::shelf_base_temp_salinity_3eqn(PetscReal gat, PetscReal sal_ocean,
//                                                        PetscReal temp_insitu, PetscReal zice, PetscReal &temp_base,
//                                                        PetscReal &sal_base){

  // The three-equation model of ice-shelf ocean interaction (Hellmer and Olbers, 1989).
  // Code derived from BRIOS subroutine iceshelf (which goes back to H.Hellmer's 2D ice shelf model code)
  // and adjusted for use in FESOM by Ralph Timmermann, 16.02.2011
  // adapted for PISM by matthias.mengel@pik-potsdam.de

  //PetscErrorCode ierr;
  PetscReal rhor, heat_flux, water_flux;
  PetscReal gats1, gats2, gat;
  // PetscReal gats1, gats2, gas, gat;
  PetscReal ep1,ep2,ep3,ep4,ep5;
  PetscReal ex1,ex2,ex3,ex4,ex5;
  PetscReal vt1,sr1,sr2,sf1,sf2,tf1,tf2,tf,sf,seta,re;
  PetscInt n, n3, nk;

  PetscReal a   = -0.0575;                // Foldvik&Kvinge (1974) [°C/psu]
  PetscReal b   =  0.0901;                // [°C]
  PetscReal c   =  7.61e-4;               // [°C/m]

  PetscReal pr  =  13.8;                  //Prandtl number      [dimensionless]
  PetscReal sc  =  2432.;                 //Schmidt number      [dimensionless]
  PetscReal ak  =  2.50e-3;               //dimensionless drag coeff.
  PetscReal sak1=  sqrt(ak);
  PetscReal un  =  1.95e-6;               //kinematic viscosity [m2/s]
  PetscReal pr1 =  pow(pr,(2./3.));       //Jenkins (1991)
  PetscReal sc1 =  pow(sc,(2./3.));
  PetscReal tob=  -20.;                   //temperature at the ice surface
  //   PetscReal rhoi=  920.;               //mean ice density
  PetscReal cpw =  4180.0;                //Barnier et al. (1995)
  PetscReal lhf =  3.33e+5;               //latent heat of fusion
  PetscReal atk =  273.15;                //0 deg C in Kelvin
  //FIXME: can use PISMs surface temp for tob?
  PetscReal cpi =  152.5+7.122*(atk+tob); //Paterson:"The Physics of Glaciers"

  PetscReal L    = 334000.;               // [J/Kg]

  // Prescribe the turbulent heat and salt transfer coeff. GAT and GAS
  //gat  = 1.00e-4;   //[m/s] RG3417

  if (gamma_T_separate_set) {
      PetscReal gat_PIG_north = gat_array[0];
      PetscReal gat_PIG_south = gat_array[1];
      PetscReal gat_TG = gat_array[2];

      PetscErrorCode ierr;

      //ierr = verbPrintf(2, grid.com,
      //              "i=%i, j=%i \n", i, j); CHKERRQ(ierr);

      if (i <= 64 && j > 82) {
          gat = gat_PIG_north;
      }

      else if (i > 64 && j > 82) {
          gat = gat_PIG_south;
      }

      else {
          gat = gat_TG;
      }
  }

  else {
      gat = gat_array[0];
  }

  //gat = gat_array[0];
  //gas  = 5.05e-7;   //[m/s] RG3417

  // Calculate
  // density in the boundary layer: rhow
  // and interface pressure pg [dbar],
  // Solve a quadratic equation for the interface salinity sb
  // to determine the melting/freezing rate seta.

  // FIXME: need to calculate water density instead of const value.
  //  call fcn_density(thetao,sal,zice,rho)
  //  rhow = density_0+rho  //was rhow= rho0+rho(i,j,N)
  //  rhow = 1028.0;
  //  rhor= rhoi/rhow;

  ep1 = cpw*gat;
  ep2 = cpi*gas;
  ep3 = lhf*gas;
  ep4 = b+c*zice;
  //  ep5 = gas/rhor;
  // negative heat flux term in the ice (due to -kappa/D)
  ex1 = a*(ep1-ep2);
  ex2 = ep1*(ep4-temp_insitu)+ep2*(tob+a*sal_ocean-ep4)-ep3;
  ex3 = sal_ocean*(ep2*(ep4-tob)+ep3);
  ex4 = ex2/ex1;
  ex5 = ex3/ex1;

  sr1 = 0.25*ex4*ex4-ex5;
  sr2 = -0.5*ex4;
  sf1 = sr2+sqrt(sr1);
  tf1 = a*sf1+ep4;
  sf2 = sr2-sqrt(sr1);
  tf2 = a*sf2+ep4;

  // Salinities < 0 psu are not defined, therefore pick the positive of the two solutions:
  if(sf1 > 0.){
    tf = tf1;
    sf = sf1;
  }else{
    tf = tf2;
    sf = sf2;
  }

  temp_base = tf;
  sal_base  = sf;

  return 0;
}

PetscErrorCode POGivenTH::compute_meltrate_3eqn( vector<double> gat_array, PetscReal gas,
						 PetscReal rhow, PetscReal rhoi,
                                                 PetscReal temp_base, PetscReal sal_base,
                                                 PetscReal sal_ocean, PetscReal &meltrate){

//PetscErrorCode POGivenTH::compute_meltrate_3eqn( PetscReal rhow, PetscReal rhoi,
//                                                 PetscReal temp_base, PetscReal sal_base,
//                                                 PetscReal sal_ocean, PetscReal &meltrate){

  // The three-equation model of ice-shelf ocean interaction (Hellmer and Olbers, 1989).
  // Code derived from BRIOS subroutine iceshelf (which goes back to H.Hellmer's 2D ice shelf model code)
  // and adjusted for use in FESOM by Ralph Timmermann, 16.02.2011
  // adapted for PISM by matthias.mengel@pik-potsdam.de

  //PetscErrorCode ierr;
  PetscReal rhor, heat_flux, water_flux;
  PetscReal gat;
  //PetscReal gas, gat;
  PetscReal ep5;

  PetscReal tob=  -20.;                        //temperature at the ice surface
  //   PetscReal rhoi=  920.;                      //mean ice density
  PetscReal cpw =  4180.0;                    //Barnier et al. (1995)
  PetscReal lhf =  3.33e+5;                   //latent heat of fusion
  PetscReal atk =  273.15;                    //0 deg C in Kelvin
  //FIXME: can use PISMs surface temp for tob?
  PetscReal cpi =  152.5+7.122*(atk+tob);     //Paterson:"The Physics of Glaciers"

  PetscReal L    = 334000.;                   // [J/Kg]

  // Prescribe the turbulent heat and salt transfer coeff. GAT and GAS

  //gat  = 1.00e-4;   //[m/s] RG3417
  //if (gamma_T_separate_set) {
  //    PetscReal gat_PIG_north = gat_array[0];
  //    PetscReal gat_PIG_south = gat_array[1];
  //    PetscReal gat_TG = gat_array[2];
  //}

  //gas  = 5.05e-7;   //[m/s] RG3417

  // Calculate
  // density in the boundary layer: rhow
  // and interface pressure pg [dbar],
  // Solve a quadratic equation for the interface salinity sb
  // to determine the melting/freezing rate seta.

  // FIXME: need to calculate water density instead of const value.
  //  call fcn_density(thetao,sal,zice,rho)
  //  rhow = density_0+rho  //was rhow= rho0+rho(i,j,N)
  //  rhow = 1028.0;
  rhor= rhoi/rhow;
  ep5 = gas/rhor;

  // Calculate the melting/freezing rate [m/s]
  meltrate = ep5*(1.0-sal_ocean/sal_base);     //rt thinks this is not needed

  //   heat_flux  = rhow*cpw*gat*(temp_insitu-tf);      // [W/m2]
  //   water_flux =          gas*(sf-sal)/sf;   // [m/s]

  return 0;
}

//FIXME: Is this method really needed???
/* This method is called by pttmpr(), which is called by potit(), which is called
   by calculate_boundlayer_temp_and_salt().*/
PetscErrorCode POGivenTH::adlprt(PetscReal salz,PetscReal temp_insitu, PetscReal pres, PetscReal &adlprt_out){
  // Berechnet aus dem Salzgehalt/psu (SALZ), der in-situ Temperatur/degC
  // (TEMP) und dem in-situ Druck/dbar (PRES) den adiabatischen Temperatur-
  // gradienten/(K Dbar^-1) ADLPRT.
  // Checkwert: ADLPRT =     3.255976E-4 K dbar^-1
  //       fuer SALZ   =    40.0 psu
  //            TEMP   =    40.0 DegC
  //            PRES   = 10000.000 dbar

  PetscReal ds;
  const PetscReal s0 = 35.0;
  const PetscReal a0 = 3.5803e-5,   a1 = 8.5258e-6,   a2 = -6.8360e-8, a3 = 6.6228e-10;
  const PetscReal b0 = 1.8932e-6,   b1 = -4.2393e-8;
  const PetscReal c0 = 1.8741e-8,   c1 = -6.7795e-10, c2 = 8.7330e-12, c3 = -5.4481e-14;
  const PetscReal d0 = -1.1351e-10, d1 = 2.7759e-12;
  const PetscReal e0 = -4.6206e-13, e1 = 1.8676e-14,  e2 = -2.1687e-16;

  ds = salz-s0;
  adlprt_out = (( ( (e2*temp_insitu + e1)*temp_insitu + e0 )*pres + ( (d1*temp_insitu + d0)*ds
                                                                      + ( (c3*temp_insitu + c2)*temp_insitu + c1 )*temp_insitu + c0 ) )*pres
                + (b1*temp_insitu + b0)*ds +  ( (a3*temp_insitu + a2)*temp_insitu + a1 )*temp_insitu + a0);

  return 0;
}

PetscErrorCode POGivenTH::pttmpr(PetscReal salz,PetscReal temp_insitu,PetscReal pres,PetscReal rfpres,
                                 PetscReal& thetao){
  // Berechnet aus dem Salzgehalt/psu (SALZ), der in-situ Temperatur/degC
  // (TEMP) und dem in-situ Druck/dbar (PRES) die potentielle Temperatur/
  // degC (PTTMPR) bezogen auf den Referenzdruck/dbar (RFPRES). Es wird
  // ein Runge-Kutta Verfahren vierter Ordnung verwendet.
  // Checkwert: PTTMPR = 36.89073 DegC
  //       fuer SALZ   =    40.0 psu
  //            TEMP   =    40.0 DegC
  //            PRES   = 10000.000 dbar
  //            RFPRES =     0.000 dbar

  PetscReal ct2  = 0.29289322 , ct3  = 1.707106781;
  PetscReal cq2a = 0.58578644 , cq2b = 0.121320344;
  PetscReal cq3a = 3.414213562, cq3b = -4.121320344;

  //real salz,temp_insitu,pres,rfpres
  PetscReal p,t,dp,dt,q, dd;
  //real adlprt

  p  = pres;
  t  = temp_insitu;
  dp = rfpres-pres;
  adlprt(salz,t,p,dd);
  dt = dp*dd;
  t  = t +0.5*dt;
  q = dt;
  p  = p +0.5*dp;
  adlprt(salz,t,p,dd);
  dt = dp*dd;
  t  = t + ct2*(dt-q);
  q  = cq2a*dt + cq2b*q;
  adlprt(salz,t,p,dd);
  dt = dp*dd;
  t  = t + ct3*(dt-q);
  q  = cq3a*dt + cq3b*q;
  p  = rfpres;
  adlprt(salz,t,p,dd);
  dt = dp*dd;
  thetao = t+ (dt-q-q)/6.0;

  return 0;
}

PetscErrorCode POGivenTH::potit(PetscReal salz,PetscReal thetao,PetscReal pres,PetscReal rfpres, PetscReal &temp_insitu_out){
  // *********************************************************************
  // Berechnet aus dem Salzgehalt[psu] (SALZ), der pot. Temperatur[oC]
  // (PT) und dem Referenzdruck[dbar] (REFPRES) die in-situ Temperatur
  // [oC] (TIN) bezogen auf den in-situ Druck[dbar] (PRES) mit Hilfe
  // eines Iterationsverfahrens aus.

  //integer iter
  //real salz,thetao,pres,rfpres,tin
  //real epsi,tpmd,pt1,ptd,pttmpr

  //tpmd / 0.001 /
  PetscReal tpmd = 0.001, epsi = 0., tin, pt1, ptd;

  for (PetscInt iter=0; iter<101; ++iter){
    //do iter=1,100
    tin  = thetao+epsi;
    //     print "tin=" + str(tin)
    pttmpr(salz,tin,pres,rfpres,pt1);
    //     PetscErrorCode ierr = verbPrintf(2, grid.com, "pt1=%e\n", pt1); CHKERRQ(ierr);
    //     print "pt1=" + str(pt1)
    ptd  = pt1-thetao;
    //     print "ptd=" + str(ptd)
    if(PetscAbs(ptd) < tpmd){
      break;
    }else{
      epsi = epsi-ptd;
    }
    if(iter==100){ SETERRQ(grid.com, 1, "in situ temperature calculation not converging."); }
  }

  temp_insitu_out = tin;
  return 0;
}


