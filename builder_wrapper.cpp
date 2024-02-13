//#include "PosData.hpp"
//
//extern "C" { 
//# include "rtklib.h"
//};
//
//extern PosDataBuilder posDataBuilder;
//extern "C" void buildDatPos(gtime_t t, const export_t* data);
//
//
//void buildDatPos(gtime_t t, const export_t* data) {
//
//	posDataBuilder.clearPosData();
//	PosDataBuilder::PosData& posData = posDataBuilder.setPosData();
//
//	posData.sysTime = t;
//	posData.nTotSats = data->ntot;
//
//	for (int i = 0; i < posData.nTotSats; i++) {
//
//		double* pr = data->range + i*2;
//		// Copy range
//		if (pr[0] != 0.0){
//			posData.range[0][i] = pr[0];
//			posData.rangeAvail[0][i] = true;
//		}	
//		else 
//			posData.rangeAvail[0][i] = false;
//
//		if (pr[1] != 0.0) {
//			posData.range[1][i] = pr[1];
//			posData.rangeAvail[1][i] = true;
//		}
//		else
//			posData.rangeAvail[1][i] = false;
//		
//		// Copy phase
//		double* pp = data->phase + i*2;
//		if (pp[0] != 0.0) {
//			posData.phase[0][i] = pp[0];
//			posData.phaseAvail[0][i] = true;
//		}
//		else
//			posData.phaseAvail[0][i] = false;
//
//		if (pp[1] != 0.0) {
//			posData.phase[1][i] = pp[1];
//			posData.phaseAvail[1][i] = true;
//		}
//		else
//			posData.phaseAvail[1][i] = false;
//
//		// Copy doppler
//		double* pd = data->doppler + i*2;
//		if (pd[0] != 0.0) {
//			posData.doppler[0][i] = pd[0];
//			posData.dopplerAvail[0][i] = true;
//		}
//		else
//			posData.dopplerAvail[0][i] = false;
//
//		if (pd[1] != 0.0) {
//			posData.doppler[1][i] = pd[1];
//			posData.dopplerAvail[1][i] = true;
//		}
//		else
//			posData.dopplerAvail[1][i] = false;
//
//		// Copy satpos
//		double* psp = data->satpos + i*6;
//		posData.satPos[i][0] = psp[0];
//		posData.satPos[i][1] = psp[1];
//		posData.satPos[i][2] = psp[2];
//
//		// Copy dts
//		double* pdts = data->satclk + i*2;
//		posData.satClk[i] = pdts[0];
//		posData.satClkDt[i] = pdts[1];
//
//		// Prn ans system
//		int prn = 0;
//		int sys = satsys(data->sat[i], &prn);
//		posData.prn[i] = prn;
//
//		switch (sys) {
//		case SYS_NONE:
//			posData.navsys[i] = NavSys::UNKNOWN;
//			break;
//		case SYS_GPS:
//			posData.navsys[i] = NavSys::GPS;
//			break;
//		case SYS_GLO:
//			posData.navsys[i] = NavSys::GLN;
//			break;
//		case SYS_GAL:
//			posData.navsys[i] = NavSys::GAL;
//			break;
//		case SYS_QZS:
//			posData.navsys[i] = NavSys::QZSS;
//			break;
//		case SYS_CMP:
//			posData.navsys[i] = NavSys::BDS;
//			break;
//		}
//
//		// Elevation, azimuth 
//		posData.azimuth[i] = data->azel[0];
//		posData.elev[i] = data->azel[1];
//	
//	}
//
//	// posDataBuilder.serialize2Json();
//
//}
