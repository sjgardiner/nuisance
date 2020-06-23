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

#ifndef MICROBOONE_SIGNALDEF_H_SEEN
#define MICROBOONE_SIGNALDEF_H_SEEN

#include "SignalDef.h"

class TSpline3;

namespace SignalDef {
  namespace MicroBooNE {

/**
 * CC with 1 muon + N>0 protons
 */
bool isCC1MuNp(FitEvent* event, double EnuMin, double EnuMax);

/**
 * CCQE-like, arxiv:2006.00108
 */
bool isCCQE(FitEvent* event, double EnuMin, double EnuMax, bool fullPS, TSpline3* prange);

  }  // namespace MicroBooNE
}  // namespace SignalDef

#endif

