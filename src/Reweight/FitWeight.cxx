#include "FitWeight.h"

#include "GENIEWeightEngine.h"
#include "LikelihoodWeightEngine.h"
#include "ModeNormEngine.h"
#include "NEUTWeightEngine.h"
#include "NIWGWeightEngine.h"
#include "NUISANCEWeightEngine.h"
#include "NuWroWeightEngine.h"
#include "OscWeightEngine.h"
#include "SampleNormEngine.h"
#include "SplineWeightEngine.h"
#include "T2KWeightEngine.h"

#ifdef __NOVA_ENABLED__
#include "NOvARwgtEngine.h"
#endif

void FitWeight::AddRWEngine(int type) {
  LOG(FIT) << "Adding reweight engine " << type << std::endl;
  switch (type) {
  case kNEUT:
    fAllRW[type] = new NEUTWeightEngine("neutrw");
    break;

  case kNUWRO:
    fAllRW[type] = new NuWroWeightEngine("nuwrorw");
    break;

  case kGENIE:
    fAllRW[type] = new GENIEWeightEngine("genierw");
    break;

  case kNORM:
    fAllRW[type] = new SampleNormEngine("normrw");
    break;

  case kLIKEWEIGHT:
    fAllRW[type] = new LikelihoodWeightEngine("likerw");
    break;

  case kT2K:
    fAllRW[type] = new T2KWeightEngine("t2krw");
    break;

  case kCUSTOM:
    fAllRW[type] = new NUISANCEWeightEngine("nuisrw");
    break;

  case kSPLINEPARAMETER:
    fAllRW[type] = new SplineWeightEngine("splinerw");
    break;

  case kNIWG:
    fAllRW[type] = new NIWGWeightEngine("niwgrw");
    break;
  case kOSCILLATION:
    fAllRW[type] = new OscWeightEngine();
    break;
  case kMODENORM:
    fAllRW[type] = new ModeNormEngine();
    break;
#ifdef __NOVA_ENABLED__
  case kNOvARWGT:
    fAllRW[type] = new NOvARwgtEngine();
    break;
#endif
  default:
    THROW("CANNOT ADD RW Engine for unknown dial type: " << type);
    break;
  }
}

WeightEngineBase *FitWeight::GetRWEngine(int type) {
  if (HasRWEngine(type)) {
    return fAllRW[type];
  }
  THROW("CANNOT get RW Engine for dial type: " << type);
}

bool FitWeight::HasRWEngine(int type) {
  switch (type) {
  case kNEUT:
  case kNUWRO:
  case kGENIE:
  case kNORM:
  case kLIKEWEIGHT:
  case kT2K:
  case kCUSTOM:
  case kSPLINEPARAMETER:
  case kNIWG:
  case kOSCILLATION:
#ifdef __NOVA_ENABLED__
  case kNOvARWGT: {
    return fAllRW.count(type);
  }
#endif
  default: { THROW("CANNOT get RW Engine for dial type: " << type); }
  }
}

void FitWeight::IncludeDial(std::string name, std::string type, double val) {
  // Should register the dial here.
  int typeenum = Reweight::ConvDialType(type);
  IncludeDial(name, typeenum, val);
}

void FitWeight::IncludeDial(std::string name, int dialtype, double val) {
  // Get the dial enum
  int nuisenum = Reweight::ConvDial(name, dialtype);

  if (nuisenum == -1) {
    ERR(FTL) << "Could not include dial " << name << std::endl;
    ERR(FTL) << "With dialtype: " << dialtype << std::endl;
    ERR(FTL) << "With value: " << val << std::endl;
    THROW("With nuisenum: " << nuisenum);
  }

  // Setup RW Engine Pointer
  if (fAllRW.find(dialtype) == fAllRW.end()) {
    AddRWEngine(dialtype);
  }
  WeightEngineBase *rw = fAllRW[dialtype];

  // Include the dial
  rw->IncludeDial(name, val);

  // Set Dial Value
  if (val != -9999.9) {
    rw->SetDialValue(name, val);
  }

  // Sort Maps
  fAllEnums[name] = nuisenum;
  fAllValues[nuisenum] = val;

  // Sort Lists
  fNameList.push_back(name);
  fEnumList.push_back(nuisenum);
  fValueList.push_back(val);
}

void FitWeight::Reconfigure(bool silent) {
  // Reconfigure all added RW engines
  for (std::map<int, WeightEngineBase *>::iterator iter = fAllRW.begin();
       iter != fAllRW.end(); iter++) {
    (*iter).second->Reconfigure(silent);
  }
}

void FitWeight::SetDialValue(std::string name, double val) {
  // Add extra check, if name not found look for one with name in it.
  int nuisenum = fAllEnums[name];
  SetDialValue(nuisenum, val);
}

// Allow for name aswell using GlobalList to determine sample name.
void FitWeight::SetDialValue(int nuisenum, double val) {
  // Conv dial type
  int dialtype = Reweight::GetDialType(nuisenum);

  if (fAllRW.find(dialtype) == fAllRW.end()) {
    ERR(FTL) << "Can't find RW engine for parameter " << fNameList[dialtype] << std::endl;
    ERR(FTL) << "With dialtype " << dialtype << ", " << Reweight::RemoveDialType(nuisenum) << std::endl;
    THROW("Are you sure you enabled the right engines?");
  }

  // Get RW Engine for this dial
  fAllRW[dialtype]->SetDialValue(nuisenum, val);
  fAllValues[nuisenum] = val;

  // Update ValueList
  for (size_t i = 0; i < fEnumList.size(); i++) {
    if (fEnumList[i] == nuisenum) {
      fValueList[i] = val;
    }
  }
}

void FitWeight::SetAllDials(const double *x, int n) {
  for (size_t i = 0; i < (UInt_t)n; i++) {
    int rwenum = fEnumList[i];
    SetDialValue(rwenum, x[i]);
  }
  Reconfigure();
}

double FitWeight::GetDialValue(std::string name) {
  // Add extra check, if name not found look for one with name in it.
  int nuisenum = fAllEnums[name];
  return GetDialValue(nuisenum);
}

double FitWeight::GetDialValue(int nuisenum) { return fAllValues[nuisenum]; }

int FitWeight::GetDialPos(std::string name) {
  int rwenum = fAllEnums[name];
  return GetDialPos(rwenum);
}

int FitWeight::GetDialPos(int nuisenum) {
  for (size_t i = 0; i < fEnumList.size(); i++) {
    if (fEnumList[i] == nuisenum) {
      return i;
    }
  }
  ERR(FTL) << "No Dial Found! " << std::endl;
  throw;
  return -1;
}

bool FitWeight::DialIncluded(std::string name) {
  return (fAllEnums.find(name) != fAllEnums.end());
}

bool FitWeight::DialIncluded(int rwenum) {
  return (fAllValues.find(rwenum) != fAllValues.end());
}

double FitWeight::CalcWeight(BaseFitEvt *evt) {
  double rwweight = 1.0;
  for (std::map<int, WeightEngineBase *>::iterator iter = fAllRW.begin();
       iter != fAllRW.end(); iter++) {
    double w = (*iter).second->CalcWeight(evt);
    rwweight *= w;
  }
  return rwweight;
}

void FitWeight::UpdateWeightEngine(const double *x) {
  size_t count = 0;
  for (std::vector<int>::iterator iter = fEnumList.begin();
       iter != fEnumList.end(); iter++) {
    SetDialValue((*iter), x[count]);
    count++;
  }
}

void FitWeight::GetAllDials(double *x, int n) {
  for (int i = 0; i < n; i++) {
    x[i] = GetDialValue(fEnumList[i]);
  }
}

// bool FitWeight::NeedsEventReWeight(const double* x) {
//   bool haschange = false;
//   size_t count = 0;

//   // Compare old to new and decide if RW needed.
//   for (std::vector<int>::iterator iter = fEnumList.begin();
//        iter != fEnumList.end(); iter++) {
//     int nuisenum = (*iter);
//     int type = (nuisenum / 1000) - (nuisenum % 1000);

//     // Compare old to new
//     double oldval = GetDialValue(nuisenum);
//     double newval = x[count];
//     if (oldval != newval) {
//       if (fAllRW[type]->NeedsEventReWeight()) {
//         haschange = true;
//       }
//     }

//     count++;
//   }

//   return haschange;
// }

double FitWeight::GetSampleNorm(std::string name) {
  if (name.empty())
    return 1.0;

  // Find norm dial
  if (fAllEnums.find(name + "_norm") != fAllEnums.end()) {
    if (fAllValues.find(fAllEnums[name + "_norm"]) != fAllValues.end()) {
      return fAllValues[fAllEnums[name + "_norm"]];
    } else {
      return 1.0;
    }
  } else {
    return 1.0;
  }
}

void FitWeight::Print() {
  LOG(REC) << "Fit Weight State: " << std::endl;
  for (size_t i = 0; i < fNameList.size(); i++) {
    LOG(REC) << " -> Par " << i << ". " << fNameList[i] << " " << fValueList[i]
             << std::endl;
  }
}
