#ifndef NUISANCESyst_H
#define NUISANCESyst_H
#include "GeneralUtils.h"

namespace Reweight {

enum NUISANCESyst {
	kUnkownNUISANCEDial = 0,
	kGaussianCorr_CCQE_norm,
	kGaussianCorr_CCQE_tilt,
	kGaussianCorr_CCQE_Pq0,
	kGaussianCorr_CCQE_Wq0,
	kGaussianCorr_CCQE_Pq3,
	kGaussianCorr_CCQE_Wq3,
	kGaussianCorr_2p2h_norm,
	kGaussianCorr_2p2h_tilt,
	kGaussianCorr_2p2h_Pq0,
	kGaussianCorr_2p2h_Wq0,
	kGaussianCorr_2p2h_Pq3,
	kGaussianCorr_2p2h_Wq3,
	kGaussianCorr_2p2h_PPandNN_norm,
	kGaussianCorr_2p2h_PPandNN_tilt,
	kGaussianCorr_2p2h_PPandNN_Pq0,
	kGaussianCorr_2p2h_PPandNN_Wq0,
	kGaussianCorr_2p2h_PPandNN_Pq3,
	kGaussianCorr_2p2h_PPandNN_Wq3,
	kGaussianCorr_2p2h_NP_norm,
	kGaussianCorr_2p2h_NP_tilt,
	kGaussianCorr_2p2h_NP_Pq0,
	kGaussianCorr_2p2h_NP_Wq0,
	kGaussianCorr_2p2h_NP_Pq3,
	kGaussianCorr_2p2h_NP_Wq3,

	kNUISANCEDial_LAST
};

int ConvertNUISANCEDial(std::string type);
std::string ConvNUISANCEDial(int type);

};
#endif