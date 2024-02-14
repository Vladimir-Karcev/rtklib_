#include "rtklib.h"
#include <string>
#include <stdexcept>
#include <iostream>
#include <iomanip> 
#include <fstream>

#include <iostream>
#include <chrono>
#include <random>
#include <algorithm>


void processObsdata(obs_t* obs) {

	const auto total_obs_number = obs->n;

	if (!total_obs_number)
		return;
	sortobs(obs);
};

void save_epoch(const obsd_t* obs, std::ofstream& stream, int n, double time) {

	for (int i = 0; i < n && i < 2 * MAXOBS; i++) {
		int sys, prn;
		sys = satsys(obs[i].sat, &prn);
		std::string system = sys == SYS_GPS ? "GPS" : sys == SYS_GLO ? "GLO" : sys == SYS_GAL ? "GAL" : sys == SYS_CMP ? "CMP" : sys == SYS_QZS ? "QZS" : "UNDEF";
		stream << std::fixed << std::setprecision(5) << time << ";" << system << ";" << prn << ";" << obs[i].P[0] << ";" << obs[i].L[0] << ";" << obs[i].D[0] << std::endl;
	}
};

static void soltocov(const sol_t* sol, double* P)
{
	P[0] = sol->qr[0]; /* xx or ee */
	P[4] = sol->qr[1]; /* yy or nn */
	P[8] = sol->qr[2]; /* zz or uu */
	P[1] = P[3] = sol->qr[3]; /* xy or en */
	P[5] = P[7] = sol->qr[4]; /* yz or nu */
	P[2] = P[6] = sol->qr[5]; /* zx or ue */
}

void print_solution(std::string sol_type, const double week, const double gpst, std::ofstream& out, const rtk_t& rtkData) {

	double pos_xyz[3] = {}, pos_llh[3] = {}, qov_enu[9] = {}, vel_enu[3] = {}, P[9];
	memcpy(pos_xyz, &rtkData.sol.rr, sizeof(pos_xyz));
	ecef2pos(pos_xyz, pos_llh);
	ecef2enu(pos_llh, &rtkData.sol.rr[3], vel_enu);
	soltocov(&rtkData.sol, P);
	covenu(pos_llh, P, qov_enu);

	out  << sol_type << '\t' 
		<< week << '\t' 
		<< gpst << '\t' 
		<< std::fixed << std::setprecision(6) << rtkData.sol.rr[0] << '\t' << rtkData.sol.rr[1] << '\t' << rtkData.sol.rr[2] << '\t'
		<< pos_llh[0] << '\t' << pos_llh[1] << '\t' << pos_llh[2] << '\t'
		<< vel_enu[0] << '\t' << vel_enu[1] << '\t' << vel_enu[2] << '\t'
		<< qov_enu[0] << '\t' << qov_enu[1] << '\t' << qov_enu[2] << '\t'
		<< static_cast<unsigned int>(rtkData.sol.ns) << std::endl;
};


// Output files
extern "C" FILE * spp_out_file;
extern "C" FILE * ppp_out_file;

static pcvs_t pcvss = { 0 };        /* receiver antenna parameters */
static pcvs_t pcvsr = { 0 };        /* satellite antenna parameters */


// Precise coordinates
double korea[3] = { -3053073.371051379, 4046489.4052480957, 3858236.0464178817 };
double algo[3] = { 0.918129094625652E+06, -.434607132958640E+07, 0.456197791639255E+07 };
double graz[3] = { 0.419442352541620E+07, 0.116270299921731E+07, 0.464724559605996E+07 };
double mro[3] = { -.255663017001537E+07, 0.509713829631369E+07, -.284838464514085E+07 };

// unpack method for retrieving data in network byte,
//   big endian, order (MSB first)
// increments index i by the number of bytes unpacked
// usage:
//   int i = 0;
//   float x = unpackFloat(&buffer[i], &i);
//   float y = unpackFloat(&buffer[i], &i);
//   float z = unpackFloat(&buffer[i], &i);
float unpackFloat(const void* buf, int* i) {
	const unsigned char* b = (const unsigned char*)buf;
	uint32_t temp = 0;
	*i += 4;
	temp = ((b[3] << 24) |
		(b[2] << 16) |
		(b[1] << 8) |
		b[0]);
	return *((float*)&temp);
}


int main(int argc, char* argv[]) {

	unsigned char buffer[] = { 0x34, 0xA8, 0x1A, 0x41};

	float t;
	int ptr = 0;
	t = unpackFloat(buffer, &ptr);
	memcpy(&t, &buffer, sizeof(t));

	char ionFile[] = "../Data_for_RTKLib/INX/FugroExt_20240116_int30m_Merged.INX";
	char pcvFile[] = "../Data_for_RTKLib/ANT/igs14_2185.atx";
	char rinexEph[] = "../Data_for_RTKLib/EPH/BRDM00DLR_S_20240160000_01D_MN.rnx";
	char rinexObs[] = "../Data_for_RTKLib/OBS/NDF9P-2024-01-16.24O";
	char preciseEph[] = "../Data_for_RTKLib/OC/rt_sp30160.sp3";
	char preciseClk[] = "../Data_for_RTKLib/OC/rt_crnx0160.clk";
	char dcbFile[] = "../Data_for_RTKLib/DCB/CAS0MGXRAP_20240160000_01D_01D_DCB.BSX";
	
	// Log PPP residuals
	ppp_out_file = fopen("residuals_ppp.txt", "w");
	spp_out_file = fopen("residuals_spp.txt", "w");

	// Log processing results and errors
	std::ofstream outFile, errFile, sm_file;
	outFile.open("out_coord.txt");
	errFile.open("error.txt");
	//outFile << "WEEK" << "\t" << "GPST" << "\t" << "X" << "\t" << "Y" << "\t" << "Z" << "\t" << std::endl;

	obs_t obs{};
	nav_t navData{};
	gtime_t ts = { 0 }, te = { 0 };

	/*char tsstr[] = "2022 05 12 06 41 20";
	str2time(tsstr, 0, 19, &ts);
	char testr[] = "2022 05 12 06 48 00";
	str2time(testr, 0, 19, &te);*/

	static sta_t stas[MAXRCV]; /* station infomation */
	//std::string options = "-GC1W-GC2W-GL1W-GL2W";
	std::string options = "";

	// Read observation data
	int status = readrnxt(rinexObs, 0, ts, te, 0.01, options.c_str(),
		&obs, NULL, NULL);
	if (status == 0)
		throw std::runtime_error("No data available!");
	else if (status == -1)
		throw std::runtime_error("Error while reading RINEX file!");
	processObsdata(&obs);

	// Read Navigation data
	options = "";
	status = readrnxt(rinexEph, 0, ts, te, 0.01, options.c_str(), &obs, &navData, stas);


	if (strlen(preciseEph) != 0) {
		// Read precise ephemeris and clock 
		int opt = 1;
		navData.ne = 0;
		navData.nemax = 512;
		navData.peph = (peph_t*)malloc(sizeof(peph_t) * navData.nemax);
		readsp3(preciseEph, &navData, opt);
	}

	if (strlen(preciseClk) != 0) {
		// Read precise clock 
		navData.nc = 0;
		navData.ncmax = 1024;
		navData.pclk = (pclk_t*)malloc(sizeof(pclk_t) * (navData.ncmax));
		status = readrnxc(preciseClk, &navData);
	}

	if (strlen(ionFile) != 0) {
		// Read ionex 
		readtec(ionFile, &navData, 0);
	}

	uniqnav(&navData);

	obsd_t data[MAXOBS];
	sol_t sol = { {0} };
	int i, j, m, iobs;
	char msg[128];
	double epoch[6];

	traceopen("trace.txt");
	tracelevel(2);

	prcopt_t my_options = {										/* defaults processing options */
	/* mode,soltype,nf,navsys */
	PMODE_PPP_KINEMA,0,2, SYS_GPS | SYS_GAL | SYS_GLO | SYS_CMP,
	15.0 * D2R,{{0,0}},											/* elmin,snrmask */
	EPHOPT_PREC,0,0,0,											/* sateph,modear,glomodear,bdsmodear */
	5,0,10,1,													/* maxout,minlock,minfix,armaxiter */
	IONOOPT_TEC,TROPOPT_EST,0,1,								/* estion,esttrop,dynamics,tidecorr */
	3,0,0,0,0,													/* niter,codesmooth,intpref,sbascorr,sbassatsel */
	0,0,														/* rovpos,refpos */
	{100.0, 100.0, 100.0},										/* eratio[] */
	{100.0,0.05,0.05,0.0,1.0},									/* err[] */
	{50.0,0.05,0.5},											/* std[] */
	{1E-7,4E-2,1E-4,1E-1,1E-2,0.0},								/* prn[] */
	GLOICB_OFF,
	5E-12,														/* sclkstab */
	{3.0,0.9999,0.25,0.1,0.05},									/* thresar */
	0.0,0.0,0.05,												/* elmaskar,almaskhold,thresslip */
	30.0,10.0,30.0,												/* maxtdif,maxinno,maxgdop */
	{0},{0},{0},												/* baseline,ru,rb */
	{"",""},													/* anttype */
	{{0}},{{0}},{0}												/* antdel,pcv,exsats */
	};

	// Set SNR mask
	my_options.snrmask.ena[0] = 1;
	for (int i = 0; i < NFREQ; i++) {
		for (j = 0; j < 9; j++)
			my_options.snrmask.mask[i][j] = 30.0;
	}

	// Enable additinal options
	my_options.posopt[0] = 1; // PCV correction enable
	my_options.posopt[2] = 2;
	my_options.posopt[3] = 1; // test eclipse
	my_options.posopt[4] = 0; // RAIM enable

	//my_options.exsats[14] = 1;

	/* Read PCV values */
	readpcv(pcvFile, &pcvss);
	setpcv(obs.data[0].time, &my_options, &navData, &pcvss, &pcvsr, stas);

	/* read DCB parameters from BIA or BSX file */
	readbiaf(dcbFile, &navData);

	/* Read DCB values P1 -> C1 */
	//readdcb(dcbFile, &navData, NULL);
	//readdcb(dcbFileP1P2, &navData, NULL);
	//readdcb_mgex(dcbFile, &navData, obs.data[0].time);

	// rtk structure to hold processing results
	bool first_call = true;
	rtk_t rtkData;
	rtkinit(&rtkData, &my_options);

	gtime_t time{ 0 };

	/* Random number to phase observation test*/
	std::random_device rd;
	std::mt19937::result_type seed = rd() ^ (
		(std::mt19937::result_type)
		std::chrono::duration_cast<std::chrono::seconds>(
			std::chrono::system_clock::now().time_since_epoch()
			).count() +
		(std::mt19937::result_type)
		std::chrono::duration_cast<std::chrono::microseconds>(
			std::chrono::high_resolution_clock::now().time_since_epoch()
			).count());

	std::mt19937 gen(seed);
	std::uniform_int_distribution<int> distrib(1, 6);

	double phase_bias[33]{ 0.0 };
	for (int i = 1; i < 33; i++) {
		phase_bias[i] = distrib(gen);
	}

	// Main Processing Loop
	for (iobs = 0; (m = nextobsf(&obs, &iobs, 0)) > 0; iobs += m) {

		for (i = j = 0; i < m && i < MAXOBS; i++) {
			data[j] = obs.data[iobs + i];
			if ((satsys(data[j].sat, NULL) & my_options.navsys) && my_options.exsats[data[j].sat - 1] != 1) {
				j++;
			}
		}
		time2epoch(data->time, epoch);
		int week = 0; double gpst = 0.0;
		gpst = time2gpst(data->time, &week); // Observation data time
		time = rtkData.sol.time; // previous epoch

		/* rover position by single point positioning */
		if (!pntpos(data, j, &navData, &rtkData.opt, &rtkData.sol, NULL, rtkData.ssat, msg)) {
			if (!rtkData.opt.dynamics) {
				outsolstat(&rtkData);
			}
			print_solution(std::string("SPP"), week, gpst, outFile, rtkData);
			errFile << week << '\t' << gpst << "\t" << std::string(msg) << std::endl;
			continue;
		}
		else {
			//print_solution(std::string("SPP"), week, gpst, outFile, rtkData);
		}

		if (time.time != 0) rtkData.tt = timediff(rtkData.sol.time, time);

		/* single point positioning */
		if (rtkData.opt.mode == PMODE_SINGLE) {
			outsolstat(&rtkData);
			print_solution(std::string("SPP"), week, gpst, outFile, rtkData);
			continue;
		}

		if (my_options.mode == PMODE_PPP_STATIC || my_options.mode == PMODE_PPP_KINEMA || my_options.mode == PMODE_PPP_CODE) {

			if (rtkData.tt == 0)
				rtkData.tt = 30.0;

			pppos(&rtkData, data, j, &navData);

			if (rtkData.sol.stat == SOLQ_PPP) {

				print_solution(std::string("PPP"), week, gpst, outFile, rtkData);
				outsolstat(&rtkData);
			}
			else {
				errFile << std::string(rtkData.errbuf) << std::endl;
				outsolstat(&rtkData);
			}
		}
	}

	outFile.close();
	errFile.close();
	fclose(ppp_out_file);
	fclose(spp_out_file);
	//traceclose();
}