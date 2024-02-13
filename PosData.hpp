#pragma once

#include "rtklib.h"
#include <iostream>
#include <cstring>
#include "gnss_data.hpp"
#include "json.hpp"

using namespace nlohmann;

constexpr unsigned int NSLOTS = 2;
constexpr unsigned int MAXSATS = 32;

class PosDataBuilder
{

public:

	struct PosData {

		gtime_t sysTime;
		int prn[MAXSATS] = {0};
		NavSys navsys[MAXSATS] = { NavSys::UNKNOWN };

		double range[NSLOTS][MAXSATS] = { 0.0 };
		double rangeSmooth[NSLOTS][MAXSATS] = { 0.0 };
		double phase[NSLOTS][MAXSATS] = { 0.0 };
		double doppler[NSLOTS][MAXSATS] = { 0.0 };
		double snr[NSLOTS][MAXSATS] = { 0.0 };

		bool rangeAvail[NSLOTS][MAXSATS] = { false };
		bool rangeSmoothAvail[NSLOTS][MAXSATS] = { false };
		bool phaseAvail[NSLOTS][MAXSATS] = { false };
		bool dopplerAvail[NSLOTS][MAXSATS] = { false };

		double rmsRange[NSLOTS][MAXSATS] = { 0.0 };
		double rmsRangeSmooth[NSLOTS][MAXSATS] = { 0.0 };
		double rmsPhase[NSLOTS][MAXSATS] = { 0.0 };
		double rmsDoppler[NSLOTS][MAXSATS] = { 0.0 };

		double satPos[MAXSATS][3] = { 0.0 };
		double satClk[MAXSATS] = { 0.0 };
		double satClkDt[MAXSATS] = { 0.0 };

		double elev[MAXSATS] = { 0.0 };
		double azimuth[MAXSATS] = { 0.0 };

		int nTotSats = 0;

	};

private:

	PosData posData;

public:

	void clearPosData() {
			
			std::fill(posData.prn, posData.prn + MAXSATS, 0);
			std::fill(posData.navsys, posData.navsys + MAXSATS, NavSys::UNKNOWN);
			std::fill(posData.range[0], posData.range[0] + NSLOTS * MAXSATS, 0);
			std::fill(posData.rangeSmooth[0], posData.rangeSmooth[0] + NSLOTS * MAXSATS, 0);
			std::fill(posData.phase[0], posData.phase[0] + NSLOTS * MAXSATS, 0);
			std::fill(posData.doppler[0], posData.doppler[0] + NSLOTS * MAXSATS, 0);
			std::fill(posData.snr[0], posData.snr[0] + NSLOTS * MAXSATS, 0);
			
			std::fill(posData.rangeAvail[0], posData.rangeAvail[0] + NSLOTS * MAXSATS, static_cast<bool>(0));
			std::fill(posData.phaseAvail[0], posData.phaseAvail[0] + NSLOTS * MAXSATS, static_cast<bool>(0));
			std::fill(posData.dopplerAvail[0], posData.dopplerAvail[0] + NSLOTS * MAXSATS, static_cast<bool>(0));
			std::fill(posData.rangeSmoothAvail[0], posData.rangeSmoothAvail[0] + NSLOTS * MAXSATS, static_cast<bool>(0));

			std::fill(posData.satPos[0], posData.satPos[0] + 3 * MAXSATS, 0);
			std::fill(posData.satClk, posData.satClk + MAXSATS, 0);
			std::fill(posData.satClkDt, posData.satClkDt + MAXSATS, 0);
			
			std::fill(posData.elev, posData.elev + MAXSATS, 0);
			std::fill(posData.azimuth, posData.azimuth + MAXSATS, 0);
	}

	PosDataBuilder() {}
	const PosData& getPosData() { return posData; } 
	PosData& setPosData() { return posData; }

	void serialize2Json()
	{	
		int gpswn = 0;
		double gpst = time2gpst(posData.sysTime, &gpswn);

		/* Json parser object*/
		json j;
		j["measurements"] = json::array();
		j["gpswn"] = gpswn;
		j["gpst"] = gpst;
		j["ntotsat"] = posData.nTotSats;

		// Iterate through all available satellites for urrent epoch
		for (int i = 0; i < posData.nTotSats; i++) {

			json jsat;
			jsat["prn"] = posData.prn[i];
			std::string sys = "";
			switch (posData.navsys[i]) {
			case NavSys::GPS:
				sys = "GPS";
				break;
			case NavSys::GAL:
				sys = "GAL";
				break;
			case NavSys::GLN:
				sys = "GLN";
				break;
			case NavSys::BDS:
				sys = "BDS";
				break;
			case NavSys::QZSS:
				sys = "QZSS";
				break;
			default:
				sys = "UNKNOWN";
			}
			jsat["navsys"] = sys;
			jsat["elevation"] = posData.elev[i];
			jsat["azimuth"] = posData.azimuth[i];

			if (posData.rangeAvail[0][i])
				jsat["rangeL1"] = posData.range[0][i];
			else
				jsat["rangeL1"] = 0.0;

			if (posData.rangeAvail[1][i])
				jsat["rangeL5"] = posData.range[1][i];
			else
				jsat["rangeL5"] = 0.0;

			if (posData.phaseAvail[0][i])
				jsat["phaseL1"] = posData.phase[0][i];
			else
				jsat["phaseL1"] = 0.0;

			if (posData.phaseAvail[1][i])
				jsat["phaseL5"] = posData.phase[1][i];
			else
				jsat["phaseL5"] = 0.0;

			if (posData.dopplerAvail[0][i])
				jsat["dopplerL1"] = posData.doppler[0][i];
			else
				jsat["dopplerL1"] = 0.0;

			if (posData.dopplerAvail[1][i])
				jsat["dopplerL5"] = posData.doppler[1][i];
			else
				jsat["dopplerL5"] = 0.0;

			jsat["snrL1"] = posData.snr[0][i];
			jsat["snrL5"] = posData.snr[1][i];

			j["measurements"].push_back(jsat);
		}

		std::cout << j.dump(4) << std::endl;

	}
	
};

