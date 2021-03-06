// Copyright (C) 2010--2011 Constantine Khroulev
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

#include "SSA.hh"

void SSA::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  dict["taud"] = new SSA_taud(this, grid, *variables);
}

SSA_taud::SSA_taud(SSA *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SSA>(m, g, my_vars) {
  
  dof = 2;
  vars.resize(dof);
  // set metadata:
  vars[0].init_2d("taud_x", grid);
  vars[1].init_2d("taud_y", grid);
  
  set_attrs("X-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 1);
}

PetscErrorCode SSA_taud::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2V *taud = new IceModelVec2V;
  ierr = taud->create(grid, "taud", false); CHKERRQ(ierr);
  ierr = taud->set_metadata(vars[0], 0); CHKERRQ(ierr); 
  ierr = taud->set_metadata(vars[1], 1); CHKERRQ(ierr); 

  ierr = taud->copy_from(model->taud); CHKERRQ(ierr);

  output = taud;
  return 0;
}

