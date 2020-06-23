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
#ifndef MICROBOONE_CCQE_1D_NU_H_SEEN
#define MICROBOONE_CCQE_1D_NU_H_SEEN

#include "Measurement1D.h"
#include "TH2Poly.h"

class TSpline3;

class MicroBooNE_CCQE_XSec_1D_nu : public Measurement1D {
public:
  /// Basic Constructor
  MicroBooNE_CCQE_XSec_1D_nu(nuiskey samplekey);

  /// Virtual Destructor
  ~MicroBooNE_CCQE_XSec_1D_nu() {};

  /// Apply numu CCQE Signal Definition
  bool isSignal(FitEvent* nvect);

  /// Fill observables
  void FillEventVariables(FitEvent* customEvent);

private:
  enum Distribution { kPmu, kPp, kCosMu, kCosP };
  Distribution fDist;
  TSpline3* fProtonRange;
  bool fFullPS;
};

#endif

