// Copyright (C) 2009-2011 Andreas Aschwanden, Ed Bueler and Constantine Khroulev
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

#ifndef __DrainageCalculator_hh
#define __DrainageCalculator_hh

class NCConfigVariable;

//! Compute the rate of drainage D(omega) for temperate ice.
class DrainageCalculator {

public:
  DrainageCalculator(const NCConfigVariable &config) {
    OM1 = config.get("drainage_target_water_frac"); // 0.01
    OM2 = 2.0 * OM1;
    OM3 = 3.0 * OM1;
    DR3 = config.get("drainage_max_rate"); // 0.05 a-1 
    DR2 = 0.1 * DR3;
  }
  virtual ~DrainageCalculator() {}

  //! Return D(omega), as in figure in [\ref AschwandenBuelerKhroulevBlatter].
  virtual PetscReal get_drainage_rate(PetscReal omega) {
    if (omega > OM1) {
      if (omega > OM2) {
        if (omega > OM3) {
          return DR3;
        } else
          return DR2 + (DR3 - DR2) * (omega - OM2) / OM1;
      } else
        return DR2 * (omega - OM1) / OM1;
    } else {
      return 0.0;
    }
  }

private:
  PetscReal OM1, OM2, OM3, DR2, DR3;
};


#endif // __DrainageCalculator_hh

