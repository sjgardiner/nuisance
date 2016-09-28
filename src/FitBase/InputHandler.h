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

#ifndef INPUT_HANDLER_H
#define INPUT_HANDLER_H

#include "TObject.h"
#include "FitEvent.h"
#include "FitWeight.h"
#include "FitParameters.h"
#include "PlotUtils.h"
#include "FitUtils.h"
#include "stdint.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>

#ifdef __GENIE_ENABLED__
#include "Conventions/Units.h"
#endif

class InputHandler {

 public:
  InputHandler(){};
  ~InputHandler(){};

  InputHandler(std::string handle, std::string infile_name);
  std::string ParseInputFile(std::string inputfile);

  void ReadBinSplineFile();
  void ReadHistogramFile();
  void ReadNeutFile();
  void ReadNuanceFile();
  void ReadGenieFile();
  void ReadGiBUUFile(bool IsNuBarDominant);
  void ReadNuWroFile();
  void ReadEventSplineFile();
  void ReadJointFile();
  FitSplineHead* GetSplineHead();

  double PredictedEventRate(double low, double high, std::string intOpt="width");
  double TotalIntegratedFlux(double low, double high, std::string intOpt="width");
  FitEvent* GetEventPointer(){ return fEvent; };
  BaseFitEvt* GetSignalPointer(){ return fSignalEvent; };

  int GetNEvents(){ return this->nEvents; };
  int GetGenEvents();

  void PrintStartInput();
  void ReadEvent(unsigned int i);

  inline TH1D* GetFluxHistogram()  { return fFluxHist;   };
  inline TH1D* GetEventHistogram() { return fEventHist;  };
  inline TH1D* GetXSecHistogram()  { return fXSecHist;   };

  std::vector<TH1*> GetFluxList(){ return this->fluxList;};
  std::vector<TH1*> GetEventList(){ return this->eventList;};
  std::vector<TH1*> GetXSecList(){ return this->xsecList;};

  int GetType(){return eventType;};
  bool CanIGoFast();
  void GetTreeEntry(const Long64_t entry);
  std::string GetInputStateString();

  int eventType;
  double GetInputWeight(const int entry=-1);

 protected:

  FitEvent* fEvent;
  BaseFitEvt* fSignalEvent;
  
  FitSplineHead* fSplineHead;

  int fMaxEvents;
  int fNEvents;
  int fCurEvt;

  TH1D* fFluxHist;
  TH1D* fEventHist;
  TH1D* fXSecHist;

  std::string fName;
  std::string fInput;
  std::string fInputFile;
  std::string fInputType;
  TFile* fInputRootFile;

  std::vector<BaseFitEvt*> fAllBaseEvents;

  std::vector<int> fJointIndexLow;
  std::vector<int> fJointIndexHigh;
  std::vector<TH1D*> fJointIndexFlux;
  std::vector<double> fJointIndexScale;
  bool fIsJointInput;

  std::vector<TH1*> fXSecList;
  std::vector<TH1*> fEventList;
  std::vector<TH1*> fFluxList;

  std::vector<TArrayD> fAllSplines;
  
  // Event Objects
#ifdef __NEUT_ENABLED__
  NeutVect *fNeutVect; //!< Pointer to NEUT Events
#endif
  
#ifdef __NUWRO_ENABLED__
  event* fNuwroEvent; //!< Pointer to NuWro Events (Set to bool if NUWRO disabled)
#endif
  
#ifdef __GENIE_ENABLED__
  GHepRecord* fGenieGHep;
  NtpMCEventRecord * fGenieNtpl; 
#endif
  
#ifdef __NUANCE_ENABLED__
  NuanceEvent* fNuanceEvt;
#endif
  
  
  TChain* tn;
  FitSplineHead* splhead;

  int maxEvents;
  int nEvents;
  int curevt_i;

  // Input Event rate flux/event histograms
  TH1D* fluxHist;  //!< Flux Histogram
  TH1D* eventHist; //!< Event Histogram
  TH1D* xsecHist;  //!< XSec Histogram

  // input root files
  TFile* inRootFile; //!< Input ROOT file (e.g NEUT MC)
  std::string inFile; ///!< Name of input ROOT file
  std::string inType; ///!< Input Type

  std::vector<BaseFitEvt*> all_events;
  std::string handleName;

#ifdef __NEUT_ENABLED__
  NeutVect *neut_event; //!< Pointer to NEUT Events
#endif

#ifdef __NUWRO_ENABLED__
  event* nuwro_event; //!< Pointer to NuWro Events (Set to bool if NUWRO disabled)
#endif

#ifdef __GENIE_ENABLED__
  GHepRecord* genie_event;  //!< Pointer to GENIE GHepRecord
  NtpMCEventRecord * mcrec; //!< Pointer to GENIE NTuple Record
#endif

#ifdef __NUANCE_ENABLED__
  NuanceEvent* nuance_event;
#endif
  
  std::vector<int> joint_index_low;
  std::vector<int> joint_index_high;
  std::vector<TH1D*> joint_index_hist;
  std::vector<double> joint_index_weight;
  bool isJointInput;
  int cur_entry;

  std::vector<TH1*> xsecList;
  std::vector<TH1*> eventList;
  std::vector<TH1*> fluxList;

  std::vector<TArrayD> spline_list;
};
#endif
