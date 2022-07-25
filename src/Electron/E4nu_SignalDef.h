// Copyright 2016 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

/*******************************************************************************
 *    This file is part of NUISANCE.
 *
 *    NUISANCE is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    NUISANCE is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with NUISANCE.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

#ifndef E4NU_SIGNALDEF_H_SEEN
#define E4NU_SIGNALDEF_H_SEEN

#include "SignalDef.h"

namespace SignalDef {
  namespace e4nu {

    // Signal definition from data sets reported in Nature 599, 565–570 (2021)
    bool isEM1p0pi( FitEvent* event, double EnuMin, double EnuMax );

  }  // namespace e4nu
}  // namespace SignalDef

#endif
