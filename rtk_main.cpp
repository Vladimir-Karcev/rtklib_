//#include "rtklib.h"
//#include <string>
//#include <stdexcept>
//#include <iostream>
//#include <iomanip> 
//#include <fstream>
//
//#include <iostream>
//#include <chrono>
//#include <random>
//#include <algorithm>
//
//static int iobsu = 0;            /* current rover observation data index */
//static int iobsr = 0;            /* current reference observation data index */
//
//extern "C" FILE * spp_out_file;
//extern "C" FILE * ppp_out_file;
//
//
//static int _inputobs(obsd_t* obs, obs_t* obss, int solq, const prcopt_t* popt)
//{
//	gtime_t time = { 0 };
//	int i, nu, nr, n = 0;
//
//	trace(3, "infunc  : iobsu=%d iobsr=%d\n", iobsu, iobsr);
//
//	if (0 <= iobsu && iobsu < obss->n) {
//		settime((time = obss->data[iobsu].time));
//	}
//	/* input forward data */
//	if ((nu = nextobsf(obss, &iobsu, 1)) <= 0) return -1;
//	if (popt->intpref) {
//		for (; (nr = nextobsf(obss, &iobsr, 2)) > 0; iobsr += nr)
//			if (timediff(obss->data[iobsr].time, obss->data[iobsu].time) > -DTTOL) break;
//	}
//	else {
//		for (i = iobsr; (nr = nextobsf(obss, &i, 2)) > 0; iobsr = i, i += nr)
//			if (timediff(obss->data[i].time, obss->data[iobsu].time) > DTTOL) break;
//	}
//	nr = nextobsf(obss, &iobsr, 2);
//	if (nr <= 0) {
//		nr = nextobsf(obss, &iobsr, 2);
//	}
//	for (i = 0; i < nu && n < MAXOBS * 2; i++) obs[n++] = obss->data[iobsu + i];
//	for (i = 0; i < nr && n < MAXOBS * 2; i++) obs[n++] = obss->data[iobsr + i];
//	iobsu += nu;
//
//	return n;
//}
//
//
//void processObsdata(obs_t* obs) {
//
//	const auto total_obs_number = obs->n;
//
//	if (!total_obs_number)
//		return;
//	sortobs(obs);
//};
//
//
//static void soltocov(const sol_t* sol, double* P)
//{
//	P[0] = sol->qr[0]; /* xx or ee */
//	P[4] = sol->qr[1]; /* yy or nn */
//	P[8] = sol->qr[2]; /* zz or uu */
//	P[1] = P[3] = sol->qr[3]; /* xy or en */
//	P[5] = P[7] = sol->qr[4]; /* yz or nu */
//	P[2] = P[6] = sol->qr[5]; /* zx or ue */
//}
//
//void print_solution(std::string sol_type, const double week, const double gpst, std::ofstream& out, const rtk_t& rtkData) {
//
//	double pos_xyz[3] = {}, pos_llh[3] = {}, qov_enu[9] = {}, P[9];
//	memcpy(pos_xyz, &rtkData.sol.rr, sizeof(pos_xyz));
//	ecef2pos(pos_xyz, pos_llh);
//	soltocov(&rtkData.sol, P);
//	covenu(pos_llh, P, qov_enu);
//
//	out << std::fixed << std::setprecision(6) << sol_type << '\t' << week << '\t' << gpst << "\t"
//		<< rtkData.sol.rr[0] << "\t" << rtkData.sol.rr[1] << "\t" << rtkData.sol.rr[2] << "\t"
//		<< qov_enu[0] << "\t" << qov_enu[1] << "\t" << qov_enu[2] << "\t"
//		<< static_cast<unsigned int>(rtkData.sol.ns) << std::endl;
//};
//
//
//static pcvs_t pcvss = { 0 };        /* receiver antenna parameters */
//static pcvs_t pcvsr = { 0 };        /* satellite antenna parameters */
//
//// Precise coordinates
//double korea[3] = { -3053073.371051379, 4046489.4052480957, 3858236.0464178817 };
//double algo[3] = { 0.918129094625652E+06, -.434607132958640E+07, 0.456197791639255E+07 };
//double graz[3] = { 0.419442352541620E+07, 0.116270299921731E+07, 0.464724559605996E+07 };
//
//double base[3] = { 2864549.8393, 2186008.3217, 5245371.0733 };
//
//
//int main(int argc, char* argv[]) {
//
//	char pcv[] = "C:/Users/artem.novichkov/Desktop/work/RTK/igs14.atx";
//	char rovobs[] = "C:/Users/artem.novichkov/Desktop/work/RTK/2198.obs";
//	char baseobs[] = "C:/Users/artem.novichkov/Desktop/work/RTK/2199.obs";
//	char eph[] = "C:/Users/artem.novichkov/Desktop/work/RTK/BRDC00IGS_R_20223550000_01D_MN.rnx";
//	spp_out_file = fopen("residuals_spp.txt", "w");
//
//	// Log processing results and errors
//	std::ofstream outf, errf;
//	outf.open("out_coord.txt");
//	errf.open("error.txt");
//	traceopen("trace.txt");
//	tracelevel(2);
//
//	static int nepoch = 0;            /* number of observation epochs */
//	static int iobsu = 0;            /* current rover observation data index */
//	static int iobsr = 0;            /* current reference observation data index */
//
//	obs_t obs{};
//	nav_t navdata{};
//	gtime_t ts = { 0 };
//	gtime_t te = { 0 };
//
//	char tsstr[] = "2022 12 21 15 00 00";
//	str2time(tsstr, 0, 19, &ts);
//	char testr[] = "2022 12 21 16 30 00";
//	str2time(testr, 0, 19, &te);
//
//	static sta_t stas[MAXRCV]; /* station infomation */
//	std::string options = "";
//
//	// Read rover observation data
//	int status = readrnxt(rovobs, 1, ts, te, 0.01, options.c_str(), &obs, NULL, NULL);
//	if (status == 0)
//		throw std::runtime_error("No data available!");
//	else if (status == -1)
//		throw std::runtime_error("Error while reading RINEX file!");
//
//	// Read base observation data
//	status = readrnxt(baseobs, 2, ts, te, 0.01, options.c_str(), &obs, NULL, NULL);
//	if (status == 0)
//		throw std::runtime_error("No data available!");
//	else if (status == -1)
//		throw std::runtime_error("Error while reading RINEX file!");
//
//	processObsdata(&obs);
//
//	// Read Navigation data
//	options = "";
//	status = readrnxt(eph, 1, ts, te, 0.01, options.c_str(), &obs, &navdata, stas);
//	uniqnav(&navdata);
//
//	obsd_t data[MAXOBS];
//	sol_t sol = { {0} };
//	char msg[128];
//
//	traceopen("trace.txt");
//	tracelevel(2);
//
//	prcopt_t opts = {										   /* defaults processing options */
//	/* mode,soltype,nf,navsys */
//	PMODE_STATIC,0,2, SYS_GLO | SYS_GPS,
//	15.0 * D2R,{{0,0}},											/* elmin,snrmask */
//	EPHOPT_BRDC,1,1,1,											/* sateph,modear,glomodear,bdsmodear */
//	5,0,10,1,													/* maxout,minlock,minfix,armaxiter */
//	IONOOPT_BRDC,TROPOPT_SAAS,0,1,								/* estion,esttrop,dynamics,tidecorr */
//	3,0,0,0,0,													/* niter,codesmooth,intpref,sbascorr,sbassatsel */
//	0,0,														/* rovpos,refpos */
//	{100.0, 100.0, 100.0},										/* eratio[] */
//	{100.0,0.05,0.05,0.0,1.0},									/* err[] */
//	{50.0,0.05,0.5},											/* std[] */
//	{1E-7,4E-2,1E-4,1E-1,1E-2,0.0},								/* prn[] */
//	GLOICB_OFF,
//	5E-12,														/* sclkstab */
//	{3.0,0.9999,0.25,0.1,0.05},									/* thresar */
//	0.0,0.0,0.05,												/* elmaskar,almaskhold,thresslip */
//	30.0,10.0,30.0,												/* maxtdif,maxinno,maxgdop */
//	{0},{0},{0},												/* baseline,ru,rb */
//	{"",""},													/* anttype */
//	{{0}},{{0}},{0}												/* antdel,pcv,exsats */
//	};
//
//	// Set SNR mask
//	opts.snrmask.ena[0] = 1;
//	for (int i = 0; i < NFREQ; i++) {
//		for (int j = 0; j < 9; j++) opts.snrmask.mask[i][j] = 20.0;
//	}
//
//	for (int i = 0; i < 3; i++) opts.rb[i] = base[i];
//	opts.modear = ARMODE_INST;
//
//	// Enable additinal options
//	opts.posopt[0] = 1; // PCV correction enable
//	opts.posopt[2] = 2;
//	opts.posopt[3] = 1; // test eclipse
//	opts.posopt[4] = 0; // RAIM enable
//
//	/* Read PCV values */
//	//readpcv(pcv, &pcvss);
//	//setpcv(obs.data[0].time, &my_options, &navdata, &pcvss, &pcvsr, stas);
//
//	// rtk structure to hold processing results
//	bool first_call = true;
//	rtk_t rtkdata;
//	rtkinit(&rtkdata, &opts);
//	gtime_t time{ 0 };
//
//	int i, nobs, n, solstatic, pri[] = { 6,1,2,3,4,5,1,6 };
//
//	while ((nobs=_inputobs(data, &obs, rtkdata.sol.stat, &opts)) >= 0) {
//
//		/* exclude satellites */
//		for (i = n = 0; i < nobs; i++) {
//			if ((satsys(data[i].sat, NULL) & opts.navsys) &&
//				opts.exsats[data[i].sat - 1] != 1) data[n++] = data[i];
//		}
//		if (n <= 0) continue;
//		if (!rtkpos(&rtkdata, data, n, &navdata)) continue;
//	}
//
//	outf.close();
//	errf.close();
//	traceclose();
//}