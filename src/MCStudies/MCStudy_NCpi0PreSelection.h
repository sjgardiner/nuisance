// Copyright 2016-2021 L. Pickering, P Stowell, R. Terri, C. Wilkinson, C. Wret

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

#ifndef MCStudy_NCpi0PreSelection_H_SEEN
#define MCStudy_NCpi0PreSelection_H_SEEN
#include "Measurement1D.h"

//********************************************************************
class MCStudy_NCpi0PreSelection : public Measurement1D {
//********************************************************************

public:

  MCStudy_NCpi0PreSelection(std::string name, std::string inputfile, FitWeight *rw, std::string type, std::string fakeDataFile);
  virtual ~MCStudy_NCpi0PreSelection() {};

  //! Grab info from event
  void FillEventVariables(FitEvent *event);

  //! Define this samples signal
  bool isSignal(FitEvent *nvect);

  //! Write Files
  void Write(std::string drawOpt);

 private:
  
  double fEventScaleFactor;
  TTree* fEventTree;
  int nlep;
  int nkplus;
  int nkaon;
  double kplusmom;
  double kaonmom;
};

#endif
