/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
*
* options : -DIERS_MODEL  use IERS tide model
*           -DOUTSTAT_AMB output ambiguity parameters to solution status
*
* references :
*    [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*    [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*        2003, November 2003
*    [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*        Space Technology Library, 2004
*    [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*        May 2009
*    [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*        Code Biases, URA
*    [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*        celestial reference frames, Geophys. Res. Let., 1997
*    [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*    [8] J.Kouba, A simplified yaw-attitude model for eclipsing GPS satellites,
*        GPS Solutions, 13:1-12, 2009
*    [9] F.Dilssner, GPS IIF-1 satellite antenna phase center and attitude
*        modeling, InsideGNSS, September, 2010
*    [10] F.Dilssner, The GLONASS-M satellite yaw-attitude model, Advances in
*        Space Research, 2010
*    [11] IGS MGEX (http://igs.org/mgex)
*
* version : $Revision:$ $Date:$
* history : 2010/07/20 1.0  new
*                           added api:
*                               tidedisp()
*           2010/12/11 1.1  enable exclusion of eclipsing satellite
*           2012/02/01 1.2  add gps-glonass h/w bias correction
*                           move windupcorr() to rtkcmn.c
*           2013/03/11 1.3  add otl and pole tides corrections
*                           involve iers model with -DIERS_MODEL
*                           change initial variances
*                           suppress acos domain error
*           2013/09/01 1.4  pole tide model by iers 2010
*                           add mode of ionosphere model off
*           2014/05/23 1.5  add output of trop gradient in solution status
*           2014/10/13 1.6  fix bug on P0(a[3]) computation in tide_oload()
*                           fix bug on m2 computation in tide_pole()
*           2015/03/19 1.7  fix bug on ionosphere correction for GLO and BDS
*           2015/05/10 1.8  add function to detect slip by MW-LC jump
*                           fix ppp solutin problem with large clock variance
*           2015/06/08 1.9  add precise satellite yaw-models
*                           cope with day-boundary problem of satellite clock
*           2015/07/31 1.10 fix bug on nan-solution without glonass nav-data
*                           pppoutsolsat() -> pppoutstat()
*           2015/11/13 1.11 add L5-receiver-dcb estimation
*                           merge post-residual validation by rnx2rtkp_test
*                           support support option opt->pppopt=-GAP_RESION=nnnn
*           2016/01/22 1.12 delete support for yaw-model bug
*                           add support for ura of ephemeris
*           2018/10/10 1.13 support api change of satexclude()
*           2020/11/30 1.14 use sat2freq() to get carrier frequency
*                           use E1-E5b for Galileo iono-free LC
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_ITER    8               /* max number of iterations */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 6              /* min satellite number for solution */
#define THRES_REJECT 2.5            /* reject threshold of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0

#define VAR_POS      SQR(60.0)       /* init variance receiver position (m^2) */
#define VAR_VEL      SQR(10.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC      SQR(10.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK      SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_ZTD      SQR(0.6)        /* init variance ztd (m^2) */
#define VAR_GRA      SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB      SQR(60.0)       /* init variance dcb (m^2) */
#define VAR_BIAS     SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO     SQR(10.0)       /* init variance iono-delay */
#define VAR_GLO_IFB  SQR(0.6)        /* variance of glonass ifb */
#define VAR_IONO_RES SQR(10.0)        /* init variance iono-delay error */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

#define EFACT_GPS_L5 1.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */
#define RES_NUM		100

/* number and index of states */
//#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->ionoopt==IONOOPT_TEC?0:(opt)->nf)
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3) /* number of position solution */
#define NC(opt)     (NSYS) /* number of clock solution */
#define ND(opt)     ((opt)->mode==PMODE_PPP_CODE?0:(opt)->nf>=2?1:0) /* number of receiver DCB */
#define NICB(opt)   ((opt)->gloicb==GLOICB_OFF?0:((opt)->gloicb==GLOICB_LNF?1:((opt)->gloicb==GLOICB_QUAD?2:((opt)->gloicb==GLOICB_1SAT?NSATGLO:13))))/* number of GLONASS ICB */
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))   /* number of tropospheric parameters */
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST||(opt)->ionoopt==IONOOPT_CONST||(opt)->ionoopt==IONOOPT_ERR?MAXSAT:0)
#define NR(opt)     (NP(opt)+NC(opt)+ND(opt)+NICB(opt)+NT(opt)+NI(opt))
#define NB(opt)     ((opt)->mode==PMODE_PPP_CODE?0:NF(opt)*MAXSAT)    /* number of phase ambiguity parameters */
#define NX(opt)     (NR(opt)+NB(opt))   /* number of estimated parameters */
#define IC(s,opt)   (NP(opt)+(s))       /* state index of clocks (s=0:gps,1:glo,2:gal,3:bds) */
#define ID(opt)     (NP(opt)+NC(opt))   /* state index of receiver DCB */
#define IICB(s,opt) (NP(opt)+NC(opt)+ND(opt)+(s)-1)   /* state index of GLONASS ICB */
#define IT(opt)     (NP(opt)+NC(opt)+ND(opt)+NICB(opt))   /* state index of tropospheric parameters */
#define II(s,opt)   (NP(opt)+NC(opt)+ND(opt)+NICB(opt)+NT(opt)+(s)-1)   /* state index of ionospheric parameters */
#define IB(s,f,opt) ((opt)->mode==PMODE_PPP_CODE?0:NR(opt)+MAXSAT*(f)+(s)-1)      /* state index of phase ambiguity parameters */

residual_t ppp_residual[RES_NUM];
FILE* ppp_out_file;

/* standard deviation of state -----------------------------------------------*/
static double STD(rtk_t* rtk, int i)
{
	if (rtk->sol.stat == SOLQ_FIX) return SQRT(rtk->Pa[i + i * rtk->nx]);
	return SQRT(rtk->P[i + i * rtk->nx]);
}

#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */

//From GLAB
/*****************************************************************************
* Name        : gravitationalDelayCorrection
* Description : Obtains the gravitational delay correction for the effect of
*               general relativity (red shift) to the GPS signal
* Parameters  :
* Name                           |Da|Unit|Description
* double  *receiverPosition       I  m    Position of the receiver
* double  *satellitePosition      I  m    Position of the satellite
* Returned value (double)         O  m    Gravitational delay correction
*****************************************************************************/
extern double gravitationalDelayCorrection(const int sys, const double* receiverPosition,
	const double* satellitePosition)
{
	double	receiverModule;
	double	satelliteModule;
	double	distance;
	double  MU = MU_GPS;

	receiverModule = sqrt(receiverPosition[0] * receiverPosition[0] + receiverPosition[1] * receiverPosition[1] +
		receiverPosition[2] * receiverPosition[2]);
	satelliteModule = sqrt(satellitePosition[0] * satellitePosition[0] + satellitePosition[1] * satellitePosition[1] +
		satellitePosition[2] * satellitePosition[2]);
	distance = sqrt((satellitePosition[0] - receiverPosition[0]) * (satellitePosition[0] - receiverPosition[0]) +
		(satellitePosition[1] - receiverPosition[1]) * (satellitePosition[1] - receiverPosition[1]) +
		(satellitePosition[2] - receiverPosition[2]) * (satellitePosition[2] - receiverPosition[2]));

	switch (sys) {
	case SYS_GPS:
		MU = MU_GPS;
		break;
	case SYS_GLO:
		MU = MU_GLO;
		break;
	case SYS_GAL:
		MU = MU_GAL;
		break;
	case SYS_CMP:
		MU = MU_CMP;
		break;
	default:
		MU = MU_GPS;
		break;
	}

	return 2.0 * MU / (CLIGHT * CLIGHT) * log((satelliteModule + receiverModule + distance) / (satelliteModule + receiverModule - distance));
}

/* write solution status for PPP ---------------------------------------------*/
extern int pppoutstat(rtk_t* rtk, char* buff)
{
	ssat_t* ssat;
	double tow, pos[3], vel[3], acc[3], * x;
	int i, j, week;
	char id[32], * p = buff;

	if (!rtk->sol.stat) return 0;

	trace(3, "pppoutstat:\n");

	tow = time2gpst(rtk->sol.time, &week);

	x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;

	/* receiver position */
	p += sprintf(p, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,
		rtk->sol.stat, x[0], x[1], x[2], STD(rtk, 0), STD(rtk, 1), STD(rtk, 2));

	/* receiver velocity and acceleration */
	if (rtk->opt.dynamics) {
		ecef2pos(rtk->sol.rr, pos);
		ecef2enu(pos, rtk->x + 3, vel);
		ecef2enu(pos, rtk->x + 6, acc);
		p += sprintf(p, "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
			"%.4f,%.5f,%.5f,%.5f\n", week, tow, rtk->sol.stat, vel[0], vel[1],
			vel[2], acc[0], acc[1], acc[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	}
	/* receiver clocks */
	i = IC(0, &rtk->opt);
	p += sprintf(p, "$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
		week, tow, rtk->sol.stat, 1, x[i] * 1E9 / CLIGHT, x[i + 1] * 1E9 / CLIGHT,
		STD(rtk, i) * 1E9 / CLIGHT, STD(rtk, i + 1) * 1E9 / CLIGHT);

	/* tropospheric parameters */
	if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
		i = IT(&rtk->opt);
		p += sprintf(p, "$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow, rtk->sol.stat,
			1, x[i], STD(rtk, i));
	}
	if (rtk->opt.tropopt == TROPOPT_ESTG) {
		i = IT(&rtk->opt);
		p += sprintf(p, "$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n", week, tow,
			rtk->sol.stat, 1, x[i + 1], x[i + 2], STD(rtk, i + 1), STD(rtk, i + 2));
	}
	/* ionosphere parameters */
	if (rtk->opt.ionoopt == IONOOPT_EST) {
		for (i = 0; i < MAXSAT; i++) {
			ssat = rtk->ssat + i;
			if (!ssat->vs) continue;
			j = II(i + 1, &rtk->opt);
			if (rtk->x[j] == 0.0) continue;
			satno2id(i + 1, id);
			p += sprintf(p, "$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n", week, tow,
				rtk->sol.stat, id, rtk->ssat[i].azel[0] * R2D,
				rtk->ssat[i].azel[1] * R2D, x[j], STD(rtk, j));
		}
	}
#ifdef OUTSTAT_AMB
	/* ambiguity parameters */
	for (i = 0; i < MAXSAT; i++) for (j = 0; j < NF(&rtk->opt); j++) {
		k = IB(i + 1, j, &rtk->opt);
		if (rtk->x[k] == 0.0) continue;
		satno2id(i + 1, id);
		p += sprintf(p, "$AMB,%d,%.3f,%d,%s,%d,%.4f,%.4f\n", week, tow,
			rtk->sol.stat, id, j + 1, x[k], STD(rtk, k));
	}
#endif
	return (int)(p - buff);
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t* obs, int n, const nav_t* nav, double* rs)
{
	double rsun[3], esun[3], r, ang, erpv[5] = { 0 }, cosa;
	int i, j;
	const char* type;

	trace(3, "testeclipse:\n");

	/* unit vector of sun direction (ecef) */
	sunmoonpos(gpst2utc(obs[0].time), erpv, rsun, NULL, NULL);
	normv3(rsun, esun);

	for (i = 0; i < n; i++) {
		type = nav->pcvs[obs[i].sat - 1].type;

		if ((r = norm(rs + i * 6, 3)) <= 0.0) continue;

		/* only block IIA */
		if (*type && !strstr(type, "BLOCK IIA")) continue;

		/* sun-earth-satellite angle */
		cosa = dot(rs + i * 6, esun, 3) / r;
		cosa = cosa < -1.0 ? -1.0 : (cosa > 1.0 ? 1.0 : cosa);
		ang = acos(cosa);

		/* test eclipse */
		if (ang<PI / 2.0 || r * sin(ang)>RE_WGS84) continue;

		trace(3, "eclipsing sat excluded %s sat=%2d\n", time_str(obs[0].time, 0),
			obs[i].sat);

		for (j = 0; j < 3; j++) rs[j + i * 6] = 0.0;
	}
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
	if (fabs(beta) < 1E-12 && fabs(mu) < 1E-12) return PI;
	return atan2(-tan(beta), sin(mu)) + PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char* type, int opt, double beta, double mu,
	double* yaw)
{
	*yaw = yaw_nominal(beta, mu);
	return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char* type, int opt,
	const double* rs, double* exs, double* eys)
{
	double rsun[3], ri[6], es[3], esun[3], n[3], p[3], en[3], ep[3], ex[3], E, beta, mu;
	double yaw, cosy, siny, erpv[5] = { 0 };
	int i;

	sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

	/* beta and orbit angle */
	matcpy(ri, rs, 6, 1);
	ri[3] -= OMGE * ri[1];
	ri[4] += OMGE * ri[0];
	cross3(ri, ri + 3, n);
	cross3(rsun, n, p);
	if (!normv3(rs, es) || !normv3(rsun, esun) || !normv3(n, en) ||
		!normv3(p, ep)) return 0;
	beta = PI / 2.0 - acos(dot(esun, en, 3));
	E = acos(dot(es, ep, 3));
	mu = PI / 2.0 + (dot(es, esun, 3) <= 0 ? -E : E);
	if (mu < -PI / 2.0) mu += 2.0 * PI;
	else if (mu >= PI / 2.0) mu -= 2.0 * PI;

	/* yaw-angle of satellite */
	if (!yaw_angle(sat, type, opt, beta, mu, &yaw)) return 0;

	/* satellite fixed x,y-vector */
	cross3(en, es, ex);
	cosy = cos(yaw);
	siny = sin(yaw);
	for (i = 0; i < 3; i++) {
		exs[i] = -siny * en[i] + cosy * ex[i];
		eys[i] = -cosy * en[i] - siny * ex[i];
	}
	return 1;
}
/* phase windup model --------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char* type, int opt,
	const double* rs, const double* rr, double* phw)
{
	double exs[3], eys[3], ek[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
	double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
	int i;

	if (opt <= 0) return 1; /* no phase windup */

	/* satellite yaw attitude model */
	if (!sat_yaw(time, sat, type, opt, rs, exs, eys)) return 0;

	/* unit vector satellite to receiver */
	for (i = 0; i < 3; i++) r[i] = rr[i] - rs[i];
	if (!normv3(r, ek)) return 0;

	/* unit vectors of receiver antenna */
	ecef2pos(rr, pos);
	xyz2enu(pos, E);
	exr[0] = E[1]; exr[1] = E[4]; exr[2] = E[7]; /* x = north */
	eyr[0] = -E[0]; eyr[1] = -E[3]; eyr[2] = -E[6]; /* y = west  */

	/* phase windup effect */
	cross3(ek, eys, eks);
	cross3(ek, eyr, ekr);
	for (i = 0; i < 3; i++) {
		ds[i] = exs[i] - ek[i] * dot(ek, exs, 3) - eks[i];
		dr[i] = exr[i] - ek[i] * dot(ek, exr, 3) + ekr[i];
	}
	cosp = dot(ds, dr, 3) / norm(ds, 3) / norm(dr, 3);
	if (cosp < -1.0) cosp = -1.0;
	else if (cosp > 1.0) cosp = 1.0;
	ph = acos(cosp) / 2.0 / PI;
	cross3(ds, dr, drs);
	if (dot(ek, drs, 3) < 0.0) ph = -ph;

	*phw = ph + floor(*phw - ph + 0.5); /* in cycle */
	return 1;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, int idx, int type,
	const prcopt_t* opt)
{
	double a = opt->err[1], b = opt->err[2];
	double c = 1.0, fact = 1.0;
	double sinel = sin(el);
	
	if (opt->mode == PMODE_PPP_CODE) {
		type = 1;
	}
	else {
		c = type ? opt->err[0] : 1.0;   /* type=0:phase, 1:code */
	}
	if (sys == SYS_GPS) {
		fact = EFACT_GPS;
		if (type) 
			c = 100.0;
	}
	else if (sys == SYS_GLO) {
		fact = EFACT_GLO;
		if (type) 
			c = 100.0;
	}
	else if (sys == SYS_CMP) {
		if (type) 
			c = 500.0;
	}
	else if (sys == SYS_GAL) {
		if (type) 
			c = 100.0;
	}
	else if (sys == SYS_QZS) {
		if (type) 
			c = 100.0;
	}

	if (opt->ionoopt == IONOOPT_IFLC) fact *= 3.0;
	double varerr = SQR(fact * c) * (SQR(a) + SQR(b / sinel));

	return varerr;
}
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t* rtk, double xi, double var, int i)
{
	int j;
	rtk->x[i] = xi;
	for (j = 0; j < rtk->nx; j++) {
		rtk->P[i + j * rtk->nx] = rtk->P[j + i * rtk->nx] = (i == j) ? var : 0.0;
	}
}
/* geometry-free phase measurement -------------------------------------------*/
static double gfmeas(const obsd_t* obs, const nav_t* nav)
{
	double freq1, freq2;

	freq1 = sat2freq(obs->sat, obs->code[0], nav);
	freq2 = sat2freq(obs->sat, obs->code[1], nav);
	if (freq1 == 0.0 || freq2 == 0.0 || obs->L[0] == 0.0 || obs->L[1] == 0.0) return 0.0;
	return (obs->L[0] / freq1 - obs->L[1] / freq2) * CLIGHT;
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
static double mwmeas(const obsd_t* obs, const nav_t* nav)
{
	double freq1, freq2;

	freq1 = sat2freq(obs->sat, obs->code[0], nav);
	freq2 = sat2freq(obs->sat, obs->code[1], nav);

	if (freq1 == 0.0 || freq2 == 0.0 || obs->L[0] == 0.0 || obs->L[1] == 0.0 ||
		obs->P[0] == 0.0 || obs->P[1] == 0.0) return 0.0;
	return (obs->L[0] - obs->L[1]) * CLIGHT / (freq1 - freq2) -
		(freq1 * obs->P[0] + freq2 * obs->P[1]) / (freq1 + freq2);
}
/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(const obsd_t* obs, const nav_t* nav, const double* azel,
	const prcopt_t* opt, const double* dantr,
	const double* dants, double phw, double* L, double* P,
	double* Lc, double* Pc)
{
	double freq[NFREQ] = { 0 }, C1, C2;
	int i, ix = 0, frq, frq2, bias_ix, sys = satsys(obs->sat, NULL);

	for (i = 0; i < opt->nf; i++) {
		L[i] = P[i] = 0.0;
		/* skip if low SNR or missing observations */
		freq[i] = sat2freq(obs->sat, obs->code[i], nav);
		if (freq[i] == 0.0 || obs->L[i] == 0.0 || obs->P[i] == 0.0) continue;
		if (testsnr(0, 0, azel[1], obs->SNR[i] * SNR_UNIT, &opt->snrmask)) continue;

		/* antenna phase center and phase windup correction */
		L[i] = obs->L[i] * CLIGHT / freq[i] - dants[i] - dantr[i] - phw * CLIGHT / freq[i];
		P[i] = obs->P[i] - dants[i] - dantr[i];

		if (opt->sateph == EPHOPT_SSRAPC || opt->sateph == EPHOPT_SSRCOM) {
			/* select SSR code correction based on code */
			if (sys == SYS_GPS)
				ix = (i == 0 ? CODE_L1W - 1 : CODE_L2W - 1);
			else if (sys == SYS_GLO)
				ix = (i == 0 ? CODE_L1P - 1 : CODE_L2P - 1);
			else if (sys == SYS_GAL)
				ix = (i == 0 ? CODE_L1X - 1 : CODE_L7X - 1);
			/* apply SSR correction */
			P[i] += (nav->ssr[obs->sat - 1].cbias[obs->code[i] - 1] - nav->ssr[obs->sat - 1].cbias[ix]);
		}
		else {   /* apply code bias corrections from file */
			if (sys == SYS_GAL && (i == 1 || i == 2)) frq = 3 - i;  /* GAL biases are L1/L5 */
			else frq = i;  /* other biases are L1/L2 */
			if (frq >= MAX_CODE_BIAS_FREQ9) continue;  /* only 2 FREQ9 per system supported in code bias table */
			bias_ix = code2bias_ix(sys, obs->code[i]); /* look up bias index in table */
			if (bias_ix > 0) {  /*  0=ref code */
				P[i] += nav->cbias[obs->sat - 1][frq][bias_ix - 1]; /* code bias */
			}
		}
	}
	/* choose FREQ9 for iono-free LC */
	*Lc = *Pc = 0.0;
	frq2 = L[1] == 0 ? 2 : 1;  /* if L[1]==0, try L[2] */
	if (freq[0] == 0.0 || freq[frq2] == 0.0) return;
	C1 = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[frq2]));
	C2 = -SQR(freq[frq2]) / (SQR(freq[0]) - SQR(freq[frq2]));

	if (L[0] != 0.0 && L[frq2] != 0.0) *Lc = C1 * L[0] + C2 * L[frq2];
	if (P[0] != 0.0 && P[frq2] != 0.0) *Pc = C1 * P[0] + C2 * P[frq2];
}


/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(rtk_t* rtk, const obsd_t* obs, int n)
{
	int i, j;

	trace(3, "detslp_ll: n=%d\n", n);

	for (i = 0; i < n && i < MAXOBS; i++) for (j = 0; j < rtk->opt.nf; j++) {
		if (obs[i].L[j] == 0.0 || !(obs[i].LLI[j] & 3)) continue;
		trace(3, "detslp_ll: slip detected sat=%2d f=%d\n", obs[i].sat, j + 1);
		rtk->ssat[obs[i].sat - 1].slip[j] = 1;
	}
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
	double g0, g1;
	int i, j;

	trace(3, "detslp_gf: n=%d\n", n);

	for (i = 0; i < n && i < MAXOBS; i++) {

		if ((g1 = gfmeas(obs + i, nav)) == 0.0) continue;

		g0 = rtk->ssat[obs[i].sat - 1].gf[0];
		rtk->ssat[obs[i].sat - 1].gf[0] = g1;

		trace(4, "detslip_gf: sat=%2d gf0=%8.3f gf1=%8.3f\n", obs[i].sat, g0, g1);

		if (g0 != 0.0 && fabs(g1 - g0) > rtk->opt.thresslip) {
			trace(3, "detslip_gf: slip detected sat=%2d gf=%8.3f->%8.3f\n",
				obs[i].sat, g0, g1);

			for (j = 0; j < rtk->opt.nf; j++) rtk->ssat[obs[i].sat - 1].slip[j] |= 1;
		}
	}
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
static void detslp_mw(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
	double w0, w1;
	int i, j;

	trace(3, "detslp_mw: n=%d\n", n);

	for (i = 0; i < n && i < MAXOBS; i++) {
		if ((w1 = mwmeas(obs + i, nav)) == 0.0) continue;

		w0 = rtk->ssat[obs[i].sat - 1].mw[0];
		rtk->ssat[obs[i].sat - 1].mw[0] = w1;

		trace(4, "detslip_mw: sat=%2d mw0=%8.3f mw1=%8.3f\n", obs[i].sat, w0, w1);

		if (w0 != 0.0 && fabs(w1 - w0) > THRES_MW_JUMP) {
			trace(3, "detslip_mw: slip detected sat=%2d mw=%8.3f->%8.3f\n",
				obs[i].sat, w0, w1);

			for (j = 0; j < rtk->opt.nf; j++) rtk->ssat[obs[i].sat - 1].slip[j] |= 1;
		}
	}
}
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(rtk_t* rtk, FILE* fP)
{
	double* F, * P, * FP, * x, * xp, pos[3], Q[9] = { 0 }, Qv[9];
	int i, j, * ix, nx;

	trace(3, "udpos_ppp:\n");

	/* fixed mode */
	if (rtk->opt.mode == PMODE_PPP_FIXED) {
		for (i = 0; i < 3; i++) initx(rtk, rtk->opt.ru[i], 1E-8, i);
		return;
	}
	/* initialize position for first epoch */
	if (norm(rtk->x, 3) <= 0.0) {
		for (i = 0; i < 3; i++)
			initx(rtk, rtk->sol.rr[i], VAR_POS, i);
		if (rtk->opt.dynamics) {
			for (i = 3; i < 6; i++) initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
			for (i = 6; i < 9; i++) initx(rtk, 1E-6, VAR_ACC, i);
		}
	}
	/* static ppp mode */
	if (rtk->opt.mode == PMODE_PPP_STATIC  || rtk->opt.mode == PMODE_PPP_CODE) {
		for (i = 0; i < 3; i++) {
			rtk->P[i * (1 + rtk->nx)] += SQR(rtk->opt.prn[5]) * fabs(rtk->tt);
		}
		if (fP != NULL) {
			matfprint(rtk->P, rtk->nx, rtk->nx, 3, 3, fP);
			fprintf(fP, "%s", "---------------POS_END-------------\r\n");
		}
		return;
	}
	/* kinmatic mode without dynamics */
	if (!rtk->opt.dynamics) {
		for (i = 0; i < 3; i++) {
			initx(rtk, rtk->sol.rr[i], VAR_POS, i);
		}
		return;
	}

	/* generate valid state index */
	ix = imat(rtk->nx, 1);
	for (i = nx = 0; i < rtk->nx; i++) {
		if (rtk->x[i] != 0.0 && rtk->P[i + i * rtk->nx] > 0.0) ix[nx++] = i;
	}
	if (nx < 9) {
		free(ix);
		return;
	}
	/* state transition of position/velocity/acceleration */
	F = eye(nx); P = mat(nx, nx); FP = mat(nx, nx); x = mat(nx, 1); xp = mat(nx, 1);

	for (i = 0; i < 6; i++) {
		F[i + (i + 3) * nx] = rtk->tt;
	}


	for (i = 0; i < 3; i++) {
		F[i + (i + 6) * nx] = SQR(rtk->tt) / 2.0;
	}

	for (i = 0; i < nx; i++) {
		x[i] = rtk->x[ix[i]];
		for (j = 0; j < nx; j++) {
			P[i + j * nx] = rtk->P[ix[i] + ix[j] * rtk->nx];
		}
	}

	/* x=F*x, P=F*P*F+Q */
	matmul("NN", nx, 1, nx, 1.0, F, x, 0.0, xp);
	matmul("NN", nx, nx, nx, 1.0, F, P, 0.0, FP);
	matmul("NT", nx, nx, nx, 1.0, FP, F, 0.0, P);

	for (i = 0; i < nx; i++) {
		rtk->x[ix[i]] = xp[i];
		for (j = 0; j < nx; j++) {
			rtk->P[ix[i] + ix[j] * rtk->nx] = P[i + j * nx];
		}
	}
	/* process noise added to only acceleration */
	Q[0] = Q[4] = SQR(rtk->opt.prn[3]) * fabs(rtk->tt);
	Q[8] = SQR(rtk->opt.prn[4]) * fabs(rtk->tt);
	ecef2pos(rtk->x, pos);
	covecef(pos, Q, Qv);
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		rtk->P[i + 6 + (j + 6) * rtk->nx] += Qv[i + j * 3];
	}

	if (fP != NULL) {
		matfprint(rtk->P, rtk->nx, rtk->nx, 3, 3, fP);
		fprintf(fP, "%s", "---------------POS_END-------------\r\n");
	}


	free(ix); free(F); free(P); free(FP); free(x); free(xp);
}

/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t* rtk, FILE* fP)
{
	prcopt_t* opt = &rtk->opt;
	double dtr, tt = fabs(rtk->tt), tdif, isb_est, isb_conv;
	int i, ii, jj, sat_r, frq, ic, isys, sys = rtk->opt.navsys;

	if (sys == SYS_GLO) {
		dtr = rtk->sol.dtr[1];
		if (fabs(dtr) < 1.0e-16) dtr = 1.0e-16;
		ic = IC(1, opt);
		initx(rtk, CLIGHT * dtr, VAR_CLK, ic);

		if (opt->gloicb == GLOICB_LNF) {   //linear function of frequency number
			ic = IICB(1, opt);
			if (rtk->x[ic] == 0.0)
				initx(rtk, 0.1, VAR_CLK, ic);
		}
		return;
	}
	else if (sys == SYS_GAL) {
		dtr = rtk->sol.dtr[2];
		ic = IC(2, opt);
		if (fabs(dtr) < 1.0e-16) dtr = 1.0e-16;
		initx(rtk, CLIGHT * dtr, VAR_CLK, ic);
		return;
	}
	else if (sys == SYS_CMP) {
		dtr = rtk->sol.dtr[3];
		ic = IC(3, opt);
		// random walk process
		if (fabs(dtr) < 1.0e-16) dtr = 1.0e-16;
		// random walk process
		initx(rtk, CLIGHT * dtr, VAR_CLK, ic);
		return;
	}
	else if (sys == SYS_QZS) {
		dtr = rtk->sol.dtr[5];
		ic = IC(4, opt);
		if (fabs(dtr) < 1.0e-16) dtr = 1.0e-16;
		// random walk process
		initx(rtk, CLIGHT * dtr, VAR_CLK, ic);
		return;
	}

	double var_clk = 60.0;
	double rndwalk = 0.001;
	
	/* initialize every epoch for clock (white noise) */
	dtr = rtk->sol.dtr[0];
	if (fabs(dtr) < 1.0e-16) dtr = 1.0e-16;
	ic = IC(0, opt);
	initx(rtk, CLIGHT * dtr, var_clk, ic);

	for (i = 1; i < NSYS; i++) {
		if (!(sys & SYS_GLO) && i == 1)
			continue;
		if (!(sys & SYS_GAL) && i == 2)
			continue;
		if (!(sys & SYS_CMP) && i == 3)
			continue;
		if (!(sys & SYS_QZS) && i == 4)
			continue;

		if (sys & SYS_GLO && i == 1)
		{
			if (opt->gloicb == GLOICB_OFF || opt->gloicb == GLOICB_LNF || opt->gloicb == GLOICB_QUAD) {  //handling of ISBs

				dtr = rtk->sol.dtr[1];
				ic = IC(1, opt);

				// random walk process
				if (rtk->x[ic] == 0.0) {
					if (fabs(dtr) < 1.0e-16)
						dtr = 1.0e-16;
					initx(rtk, CLIGHT * dtr, VAR_CLK, ic);
				}
				else {
					rtk->P[ic + ic * rtk->nx] += SQR(rndwalk) * tt;
				}
			}
			if (opt->gloicb == GLOICB_LNF) {   //linear function of frequency number
				ic = IICB(1, opt);
				if (rtk->x[ic] == 0.0)
					initx(rtk, 0.1, VAR_CLK, ic);
			}
		}
		else {
			
			if (i == 4)
				dtr = rtk->sol.dtr[5]; // QZS-GPS time offset (s)
			else
				dtr = rtk->sol.dtr[i];
			
			ic = IC(i, opt);

			// White noise process for ISB
			/*if (fabs(dtr) < 1.0e-16) dtr = 1.0e-16;
			initx(rtk, CLIGHT * dtr,  3.0, ic);*/

			// Random walk process for ISB
			if (rtk->x[ic] == 0.0) {
				if (fabs(dtr) < 1.0e-16) 
					dtr = 1.0e-16;
				initx(rtk, CLIGHT * dtr, VAR_CLK, ic);
			}
			else {
				rtk->P[ic + ic * rtk->nx] += SQR(rndwalk) * tt;
			}
		}
	}

	if (fP != NULL) {
		matfprint(rtk->P, rtk->nx, rtk->nx, 3, 3, fP);
		fprintf(fP, "%s", "---------------CLK_END-------------\r\n");
	}

}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t* rtk, FILE* fP)
{
	double pos[3], azel[] = { 0.0,PI / 2.0 }, ztd, var;
	int i = IT(&rtk->opt), j;

	trace(3, "udtrop_ppp:\n");

	if (rtk->x[i] == 0.0) {
		var = SQR(0.6);
		initx(rtk, 0.15, var, i);
		if (rtk->opt.tropopt >= TROPOPT_ESTG) {
			for (j = i + 1; j < i + 3; j++) initx(rtk, 1E-6, VAR_GRA, j);
		}
	}
	else {
		rtk->P[i + i * rtk->nx] += SQR(rtk->opt.prn[2]) * fabs(rtk->tt);
		if (rtk->opt.tropopt >= TROPOPT_ESTG) {
			for (j = i + 1; j < i + 3; j++) {
				rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[2] * 0.1) * fabs(rtk->tt);
			}
		}
	}

	if (fP != NULL) {
		matfprint(rtk->P, rtk->nx, rtk->nx, 3, 3, fP);
		fprintf(fP, "%s", "-------------END_TROPO---------------\r\n");
	}

}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udiono_ppp(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav, FILE* fP)
{
	const double* lam;
	double freq1, freq2, freq3, var = 0.0, ion = 0.0, sinel = 0.0, dion_tec = 0.0, delta = 0.0, vari_tec = 0.0, elev, pos[3], * azel;
	char* p;
	int i, j, k, gap_resion = GAP_RESION;

	trace(3, "udiono_ppp:\n");

	if ((p = strstr(rtk->opt.pppopt, "-GAP_RESION="))) {
		int ret = sscanf(p, "-GAP_RESION=%d", &gap_resion);
	}
	for (i = 0; i < MAXSAT; i++) {
		j = II(i + 1, &rtk->opt);
		if (rtk->x[j] != 0.0 && (int)rtk->ssat[i].outc[0] > gap_resion) {
			rtk->x[j] = 0.0;
		}
	}
	
	ecef2pos(rtk->sol.rr, pos);
	
	for (i = 0; i < n; i++) {
		int sat = obs[i].sat;
		j = II(obs[i].sat, &rtk->opt);
		lam = nav->lam[sat - 1];
		double dion = 0.0, prec_ion = 0.0;

		if (rtk->opt.ionoopt == IONOOPT_ERR) {

			if (rtk->x[j] == 0.0) {
				k = 1;
				double diff_ion = 0.0, k_ion_var = 0.0, k_ion = 0.0;
				if (obs[i].P[k] == 0.0 || lam[k] == 0.0)
					k = 2;
				if (obs[i].P[0] == 0.0 || obs[i].P[k] == 0.0 || lam[0] == 0.0 || lam[k] == 0.0) {
					continue;
				}
				else {
					diff_ion = (obs[i].P[0] - obs[i].P[k]) / (1.0 - SQR(lam[k] / lam[0]));
					//k_ion = ionmodel(obs[i].time, nav->ion_gps, pos, rtk->ssat[sat - 1].azel);
					if (/*iontec(obs[i].time, nav, pos, rtk->ssat[sat - 1].azel, 1, &prec_ion, &vari_tec) && */ 
						diff_ion != 0.0) {
						prec_ion = ionmodel(obs[i].time, nav->ion_gps, pos, rtk->ssat[sat - 1].azel);
						ion = diff_ion - prec_ion;
					}
					else {
						continue;
					}
					var = VAR_IONO_RES;
				}
				initx(rtk, ion, var, j);
			}
			else {
				rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[1]) * fabs(rtk->tt);
			}

		} 
		else if (rtk->opt.ionoopt == IONOOPT_EST) {

			if (rtk->x[j] == 0.0) {
				k = 1;
				if (obs[i].P[k] == 0.0 || lam[k] == 0.0)
					k = 2;
				if (obs[i].P[0] == 0.0 || obs[i].P[k] == 0.0 || lam[0] == 0.0 || lam[k] == 0.0) {
					iontec(obs[i].time, nav, pos, rtk->ssat[sat - 1].azel, 1, &ion, &vari_tec);
					var = VAR_IONO;
				}
				else {
					ion = (obs[i].P[0] - obs[i].P[k]) / (1.0 - SQR(lam[k] / lam[0]));
					var = VAR_IONO;
				}
				initx(rtk, ion, var, j);
			}
			else {
				rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[1]) * fabs(rtk->tt);
			}
		}
	}

	if (fP != NULL) {
		matfprint(rtk->P, rtk->nx, rtk->nx, 3, 3, fP);
		fprintf(fP, "%s", "-------------END_IONO---------------\r\n");
	}
}
/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_ppp(rtk_t* rtk)
{
	int i = ID(&rtk->opt);

	if (rtk->x[i] == 0.0) {
		initx(rtk, 1E-6, VAR_DCB, i);
	}
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav, FILE* fP)
{
	const double* lam;
	double L[NFREQ], P[NFREQ], Lc, Pc, bias[MAXOBS], offset = 0.0, pos[3] = { 0 };
	double freq1, freq2, ion, ionvar, dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
	int i, j, k, l, f, sat, slip[MAXOBS] = { 0 }, clk_jump = 0;

	trace(3, "udbias  : n=%d\n", n);

	/* handle day-boundary clock jump */
	if (rtk->opt.posopt[5]) {
		clk_jump = ROUND(time2gpst(obs[0].time, NULL) * 10) % 864000 == 0;
	}
	for (i = 0; i < MAXSAT; i++) for (j = 0; j < rtk->opt.nf; j++) {
		rtk->ssat[i].slip[j] = 0;
	}
	/* detect cycle slip by LLI */
 	detslp_ll(rtk, obs, n);

	/* detect cycle slip by geometry-free phase jump */
	detslp_gf(rtk, obs, n, nav);

	/* detect slip by Melbourne-Wubbena linear combination jump */
	detslp_mw(rtk, obs, n, nav);

	ecef2pos(rtk->sol.rr, pos);

	int tmp = NF(&rtk->opt);
	ecef2pos(rtk->sol.rr, pos);
	for (f = 0; f < NF(&rtk->opt); f++) {

		/* reset phase-bias if expire obs outage counter */
		for (i = 0; i < MAXSAT; i++) {
			if (++rtk->ssat[i].outc[f] > (uint32_t)rtk->opt.maxout ||
				rtk->opt.modear == ARMODE_INST || clk_jump) {
				initx(rtk, 0.0, 0.0, IB(i + 1, f, &rtk->opt));
			}
		}
		for (i = k = 0; i < n && i < MAXOBS; i++) {

			sat = obs[i].sat;
			j = IB(sat, f, &rtk->opt);
			int sys = satsys(sat, NULL);

			corr_meas(obs + i, nav, rtk->ssat[sat - 1].azel, &rtk->opt, dantr, dants,
				0.0, L, P, &Lc, &Pc);

			bias[i] = 0.0;

			if (rtk->opt.ionoopt == IONOOPT_IFLC) {
				bias[i] = Lc - Pc;
				slip[i] = rtk->ssat[sat - 1].slip[0] || rtk->ssat[sat - 1].slip[1];
			}
			else if (L[f] != 0.0 && P[f] != 0.0) {

				slip[i] = rtk->ssat[sat - 1].slip[f];
				lam = nav->lam[sat - 1];

				l = 1; // Use L2 to estimate ionosphere
				if (obs[i].P[l] == 0.0 || lam[l] == 0.0)
					l = 2; // if no L2 measurements use L5
				if (obs[i].P[0] == 0.0 || obs[i].P[l] == 0.0 || lam[0] == 0.0 || lam[l] == 0.0 || lam[f] == 0.0) {
					//continue;
					iontec(obs[i].time, nav, pos, rtk->ssat[sat - 1].azel, 1, &ion, &ionvar);
				}
				else
					ion = (obs[i].P[0] - obs[i].P[l]) / (1.0 - SQR(lam[l] / lam[0]));

				bias[i] = L[f] - P[f] + 2.0 * ion * SQR(lam[f] / lam[0]);
			}

			if (rtk->x[j] == 0.0 || slip[i] || bias[i] == 0.0) continue;
			offset += bias[i] - rtk->x[j];
			k++;
		}
		/* correct phase-code jump to ensure phase-code coherency */
		if (k >= 2 && fabs(offset / k) > 0.0005 * CLIGHT) {
			for (i = 0; i < MAXSAT; i++) {
				j = IB(i + 1, f, &rtk->opt);
				if (rtk->x[j] != 0.0) rtk->x[j] += offset / k;
			}
			trace(2, "phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
				time_str(rtk->sol.time, 0), k, offset / k / CLIGHT);
		}
		for (i = 0; i < n && i < MAXOBS; i++) {
			sat = obs[i].sat;
			j = IB(sat, f, &rtk->opt);

			rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[0]) * fabs(rtk->tt);

			if (bias[i] == 0.0 || (rtk->x[j] != 0.0 && !slip[i])) continue;

			/* reinitialize phase-bias if detecting cycle slip */
			initx(rtk, bias[i], VAR_BIAS, IB(sat, f, &rtk->opt));

			/* reset fix flags */
			for (k = 0; k < MAXSAT; k++) rtk->ambc[sat - 1].flags[k] = 0;

			trace(5, "udbias_ppp: sat=%2d bias=%.3f\n", sat, bias[i]);
		}
	}

	if (fP != NULL) {
		matfprint(rtk->P, rtk->nx, rtk->nx, 3, 3, fP);
		fprintf(fP, "%s", "-------------END_PHASE_BIAS---------------\r\n");
	}

}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
	FILE* fP = NULL;
	char fileP[] = "P.txt";

	//if (!(fP = fopen(fileP, "w"))) return 0;

	trace(3, "udstate_ppp: n=%d\n", n);

	/* temporal update of position */
	udpos_ppp(rtk, fP);

	/* temporal update of clock */
	udclk_ppp(rtk, fP);

	/* temporal update of tropospheric parameters */
	if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
		udtrop_ppp(rtk, fP);
	}
	/* temporal update of ionospheric parameters */
	if (rtk->opt.ionoopt == IONOOPT_EST || rtk->opt.ionoopt == IONOOPT_ERR) {
		udiono_ppp(rtk, obs, n, nav, fP);
	}
	/* temporal update of L5-receiver-dcb parameters */
	if (rtk->opt.nf >= 5) {
		uddcb_ppp(rtk);
	}
	if (rtk->opt.mode != PMODE_PPP_CODE) {
		/* temporal update of phase-bias */
		udbias_ppp(rtk, obs, n, nav, fP);
	}
	
	trace(15, "xp_udstate=\n"); tracemat(5, rtk->x, rtk->nx, 1, 15, 5);
	trace(15, "P_udstate=\n"); tracemat(5, rtk->P, rtk->nx, rtk->nx, 15, 5);

	//fclose(fP);
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(int sat, const double* rs, const double* rr, const pcv_t* pcv,
	double* dant)
{
	double ru[3], rz[3], eu[3], ez[3], nadir, cosa;
	int i;

	for (i = 0; i < 3; i++) {
		ru[i] = rr[i] - rs[i];
		rz[i] = -rs[i];
	}
	if (!normv3(ru, eu) || !normv3(rz, ez)) return;

	cosa = dot(eu, ez, 3);
	cosa = cosa < -1.0 ? -1.0 : (cosa > 1.0 ? 1.0 : cosa);
	nadir = acos(cosa);

	antmodel_s(pcv, nadir, dant);
	//antmodel_s(pcv,nadir,dant);
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double* pos,
	const double* azel, const double* x, double* dtdx,
	double* var)
{
	const double zazel[] = { 0.0,PI / 2.0 };
	double zhd, m_h, m_w, cotz, grad_n, grad_e;

	/* zenith hydrostatic delay */
	//zhd = tropmodel(time, pos, zazel, 0.0);

	double geohgt = geoidh(pos);
	double orthohgt = pos[2] - geohgt;
	double hzd, wzd;
	double doy = time2doy(time);

	unb3tropmodel(pos, azel, orthohgt, doy, azel[0], &zhd, &wzd);

	/* mapping function */
	m_h = tropmapf(time, pos, azel, &m_w);

	if (azel[1] > 0.0) {

		/* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
		cotz = 1.0 / tan(azel[1]);
		grad_n = m_w * cotz * cos(azel[0]);
		grad_e = m_w * cotz * sin(azel[0]);
		m_w += grad_n * x[1] + grad_e * x[2];
		dtdx[1] = grad_n * (x[0] - zhd);
		dtdx[2] = grad_e * (x[0] - zhd);
	}
	dtdx[0] = m_w;
	*var = SQR(0.01);
	//return m_h * zhd + m_w * (x[0] - zhd);
	return m_h * zhd + m_w * (x[0]);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(gtime_t time, const double* pos, const double* azel,
	const prcopt_t* opt, const double* x, double* dtdx,
	const nav_t* nav, double* dtrp, double* var)
{
	double trp[3] = { 0 };

	if (opt->tropopt == TROPOPT_SAAS) {
		*dtrp = tropmodel(time, pos, azel, REL_HUMI);
		*var = SQR(ERR_SAAS);
		return 1;
	}
	if (opt->tropopt == TROPOPT_SBAS) {
		*dtrp = sbstropcorr(time, pos, azel, var);
		return 1;
	}
	if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
		matcpy(trp, x + IT(opt), opt->tropopt == TROPOPT_EST ? 1 : 3, 1);
		*dtrp = trop_model_prec(time, pos, azel, trp, dtdx, var);
		return 1;
	}
	return 0;
}
/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double* pos, const double* azel,
	const prcopt_t* opt, int sat, const double* x,
	const nav_t* nav, double* dion, double* var)
{
	if (opt->ionoopt == IONOOPT_SBAS) {
		return sbsioncorr(time, nav, pos, azel, dion, var);
	}
	if (opt->ionoopt == IONOOPT_TEC) {
		return iontec(time, nav, pos, azel, 1, dion, var);
	}
	if (opt->ionoopt == IONOOPT_BRDC) {
		*dion = ionmodel(time, nav->ion_gps, pos, azel);
		*var = SQR(*dion * ERR_BRDCI);
		return 1;
	}
	if (opt->ionoopt == IONOOPT_EST || opt->ionoopt == IONOOPT_ERR) {
		*dion = x[II(sat, opt)];
		*var = 0.0;
		return 1;
	}
	if (opt->ionoopt == IONOOPT_IFLC) {
		*dion = *var = 0.0;
		return 1;
	}
	return 0;
}
/* phase and code residuals --------------------------------------------------*/
static int ppp_res(int post, const obsd_t* obs, int n, const double* rs,
	const double* dts, const double* var_rs, const int* svh,
	const double* dr, int* exc, const nav_t* nav,
	const double* x, rtk_t* rtk, double* v, double* H, double* R,
	double* azel)
{
	const double* lam;
	prcopt_t* opt = &rtk->opt;
	double y, r, cdtr, bias, C = 0.0, C1 = 0.0, C2 = 0.0, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], Lc, Pc;
	double var[MAXOBS * 2], dtrp = 0.0, vtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0, dcb, freq;
	double dion_tec = 0.0, vari_tec = 0.0, delta_iono = 0.0;
	double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
	double ve[MAXOBS * 2 * NFREQ] = { 0 }, vmax = 0;
	char str[32];
	int ne = 0, obsi[MAXOBS * 2 * NFREQ] = { 0 }, frqi[MAXOBS * 2 * NFREQ], maxobs, maxfrq, rej;
	int i, j, k, kk, ic = 0, id = 0, sat, sys, frq, nv = 0, nx = rtk->nx, stat = 1, res_cnt = 0 , prn = 0;
	double gravitationalDelayModel = 0.0;

	time2str(obs[0].time, str, 2);

	for (i = 0; i < MAXSAT; i++)
		for (j = 0; j < opt->nf; j++)
			rtk->ssat[i].vsat[j] = 0;

	for (i = 0; i < 3; i++) rr[i] = x[i] + dr[i];
	ecef2pos(rr, pos);

	for (i = 0; i < n && i < MAXOBS; i++) {
		sat = obs[i].sat;
		lam = nav->lam[sat - 1];
		double gpst = 0.0, week = 0.0;
		gpst = time2gpst(obs[i].time, &week);

		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 ||
			satazel(pos, e, azel + i * 2) < opt->elmin) {
			exc[i] = 1;
			continue;
		}
		int sat_excluded = satexclude(obs[i].sat, var_rs[i], svh[i], opt);
		if (!(sys = satsys(sat, &prn)) || !rtk->ssat[sat - 1].vs || sat_excluded || exc[i]) {
			exc[i] = 1;
			continue;
		}

		/* tropospheric and ionospheric model */
		int tropo_status = 0;
		
		if (rtk->opt.tropopt == UNB3M) {
			tropo_status = tropcorr(obs[i].time, nav, pos, azel + i * 2, opt->tropopt, &dtrp, &vtrp);
		}
		else {
			tropo_status = model_trop(obs[i].time, pos, azel + i * 2, opt, x, dtdx, nav, &dtrp, &vart);
		}

		int iono_status = model_iono(obs[i].time, pos, azel + i * 2, opt, sat, x, nav, &dion, &vari);

		//ionospheric delay derived from GIM
		if (rtk->opt.ionoopt == IONOOPT_ERR) {
			//iontec(obs[i].time, nav, pos, azel + i * 2, 1, &dion_tec, &vari_tec);
			dion_tec = ionmodel(obs[i].time, nav->ion_gps, pos, rtk->ssat[sat - 1].azel);
			dion_tec *= SQR(lam[0] / lam_carr[0]);
		}

		if (!tropo_status || !iono_status) {
			continue;
		}

		/* satellite and receiver antenna model */
		if (opt->posopt[0])
			satantpcv(sat, rs + i * 6, rr, nav->pcvs + sat - 1, dants);

		antmodel(opt->pcvr, opt->antdel[0], azel + i * 2, opt->posopt[1], dantr);

		if (opt->mode != PMODE_PPP_CODE) {
			/* phase windup model */
			if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type,
				opt->posopt[2] ? 2 : 0, rs + i * 6, rr, &rtk->ssat[sat - 1].phw)) {
				continue;
			}
		}
		else {
			rtk->ssat[sat - 1].phw = 0.0;
		}
		
		//Gravitational delay correction */
		gravitationalDelayModel = gravitationalDelayCorrection(sys, rr, rs + i * 6);

		/* corrected phase and code measurements */
		corr_meas(obs + i, nav, azel + i * 2, &rtk->opt, dantr, dants,
			rtk->ssat[sat - 1].phw, L, P, &Lc, &Pc);

		/* stack phase and code residuals {L1,P1,L2,P2,...} */
		int obsn = 0;
		if (opt->mode != PMODE_PPP_CODE) {
			obsn = 2 * NF(opt);
		} 
		else {
			obsn = NF(opt);
		}

		for (j = 0; j < obsn; j++) {

			dcb = bias = delta_iono = 0.0;
			if (opt->ionoopt == IONOOPT_IFLC) {
				if ((y = j % 2 == 0 ? Lc : Pc) == 0.0) continue;
			}
			else if (opt->mode == PMODE_PPP_CODE) {
				if ((y = P[j]) == 0.0) continue;
			}
			else {
				if ((y = j % 2 == 0 ? L[j / 2] : P[j / 2]) == 0.0) continue;
			}

			if (opt->mode == PMODE_PPP_CODE)
				freq = sat2freq(sat, obs[i].code[j], nav);
			else 
				freq = sat2freq(sat, obs[i].code[j / 2], nav);
			
			if (freq == 0.0)
				continue;

			if (opt->mode == PMODE_PPP_CODE)
				C = SQR(lam[j] / lam[0]) * 1.0;
			else 
				C  = SQR(lam[j / 2] / lam[0]) * (j % 2 == 0 ? -1.0 : 1.0);
			
			C1 = SQR(lam[1]) / (SQR(lam[1]) - SQR(lam[0]));
			C2 = -SQR(lam[0]) / (SQR(lam[1]) - SQR(lam[0]));

			for (k = 0; k < nx; k++) H[k + nx * nv] = k < 3 ? -e[k] : 0.0;

			/* receiver clock */
			cdtr = 0.0;
			if (sys == SYS_GPS) {
				ic = IC(0, opt);
				cdtr = x[ic];
				H[ic + nx * nv] = 1.0;
			}
			if (sys == SYS_GLO) {
				if (opt->navsys == SYS_GLO || opt->gloicb == GLOICB_OFF || opt->gloicb == GLOICB_LNF || opt->gloicb == GLOICB_QUAD) {  //handling of ISBs
					ic = IC(0, opt);
					id = IC(1, opt);
					cdtr = x[ic] + x[id];
					H[ic + nx * nv] = 1.0;
					H[id + nx * nv] = 1.0;
				}
				else if (opt->navsys != SYS_GLO) {
					ic = IC(0, opt);
					cdtr = x[ic];
					H[ic + nx * nv] = 1.0;
				}
				if (opt->gloicb == GLOICB_LNF) {   //linear function of frequency number
					if ((opt->nf == 2 && (j / 2 == 1 && j % 2 == 1)) || (opt->nf == 1 && (j / 2 == 0 && j % 2 == 1))) {
						frq = get_glo_fcn(sat, nav);
						ic = IICB(1, opt);
						cdtr += frq * x[ic];
						H[ic + nx * nv] = frq;
					}
				}
			}
			if (sys == SYS_GAL) {
				ic = IC(0, opt);
				id = IC(2, opt);
				cdtr = x[ic] + x[id];
				H[ic + nx * nv] = 1.0;
				H[id + nx * nv] = 1.0;
			}
			if (sys == SYS_CMP) {
				ic = IC(0, opt);
				id = IC(3, opt);
				cdtr = x[ic] + x[id];
				H[ic + nx * nv] = 1.0;
				H[id + nx * nv] = 1.0;
			}
			if (sys == SYS_QZS) {
				ic = IC(0, opt);
				id = IC(4, opt);
				cdtr = x[ic] + x[id];
				H[ic + nx * nv] = 1.0;
				H[id + nx * nv] = 1.0;
			}

			if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
				for (k = 0; k < (opt->tropopt >= TROPOPT_ESTG ? 3 : 1); k++) {
					H[IT(opt) + k + nx * nv] = dtdx[k];
				}
			}
			if (opt->ionoopt == IONOOPT_EST || opt->ionoopt == IONOOPT_ERR) {
				int ii = II(sat, opt);
				if (x[ii] == 0.0) continue;
				H[II(sat, opt) + nx * nv] = C;
			}

			//if (sys != SYS_GLO && (j % 2 == 1)) {
			//	dcb += rtk->x[ID(opt)];
			//	H[ID(opt) + nx * nv] = 1.0;
			//}

			//if (j / 2 == 2 && j % 2 == 1) { /* L5-receiver-dcb */
			//	dcb += rtk->x[ID(opt)];
			//	H[ID(opt) + nx * nv] = 1.0;
			//}

			if ((j % 2 == 0) && opt->mode != PMODE_PPP_CODE) { /* phase bias */
				int ib = IB(sat, j / 2, opt);
				if ((bias = x[IB(sat, j / 2, opt)]) == 0.0) continue;
				H[IB(sat, j / 2, opt) + nx * nv] = 1.0;
			}

			/* residual */
			double resid = y - (r + cdtr - CLIGHT * dts[i * 2] + dtrp + C * (dion_tec + dion) + dcb + bias);
			v[nv] = resid;

			if (opt->mode == PMODE_PPP_CODE) {
				rtk->ssat[sat - 1].resp[j] = v[nv];
			}
			else {
				if (j % 2 == 0)
					rtk->ssat[sat - 1].resc[j / 2] = v[nv];
				else
					rtk->ssat[sat - 1].resp[j / 2] = v[nv];
			}

			/* variance */
			double verr = 0.0;
			verr = varerr(obs[i].sat, sys, azel[1 + i * 2], j / 2, j % 2, opt) + vart + SQR(C) * vari + var_rs[i];

			if (opt->mode == PMODE_PPP_CODE) {
				if (sys == SYS_GLO) verr += VAR_GLO_IFB;
			}
			else {
				if (sys == SYS_GLO && j % 2 == 1) var[nv] += VAR_GLO_IFB;
			}

			var[nv] = verr;

			trace(4, "%s sat=%2d %s%d res=%9.4f sig=%9.4f el=%4.1f\n", str, sat,
				j % 2 ? "P" : "L", j / 2 + 1, v[nv], sqrt(var[nv]), azel[1 + i * 2] * R2D);

			/* reject satellite by pre-fit residuals */
			if (!post && opt->maxinno > 0.0 && fabs(v[nv]) > opt->maxinno) {
				trace(2, "outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
					post, str, sat, j % 2 ? "P" : "L", j / 2 + 1, v[nv], azel[1 + i * 2] * R2D);
				exc[i] = 1; 

				if (opt->mode == PMODE_PPP_CODE) {
					rtk->ssat[sat - 1].rejc[1]++;
				}
				else {
					rtk->ssat[sat - 1].rejc[j % 2]++;
				}
				continue;
			}
			/* record large post-fit residuals */
			if (post && fabs(v[nv]) > sqrt(var[nv]) * THRES_REJECT) {
				obsi[ne] = i; frqi[ne] = j; ve[ne] = v[nv]; ne++;
			}
			if (opt->mode == PMODE_PPP_CODE) {
				rtk->ssat[sat - 1].vsat[j] = 1;
			}
			else{
				if (j % 2 == 0) rtk->ssat[sat - 1].vsat[j / 2] = 1;
			}
			
			nv++;

			// record residual
			{
				int idx = res_cnt;
				ppp_residual[idx].sat = sat;
				ppp_residual[idx].sys = sys;
				ppp_residual[idx].prn = prn;
				ppp_residual[idx].freq = j;
				ppp_residual[idx].gpst = gpst;
				ppp_residual[idx].week = week;
				ppp_residual[idx].prange = obs[i].P[j];
				ppp_residual[idx].doppler = obs[i].D[j];
				ppp_residual[idx].phase = obs[i].L[j];
				ppp_residual[idx].az = azel[i * 2] * R2D;
				ppp_residual[idx].el = azel[1 + i * 2] * R2D;
				ppp_residual[idx].snr = obs[i].SNR[j] * SNR_UNIT;
				ppp_residual[idx].range = r;
				ppp_residual[idx].dts = CLIGHT * dts[i * 2];
				ppp_residual[idx].dtr = cdtr;
				ppp_residual[idx].glo_isb = x[IC(1, opt)];
				ppp_residual[idx].gal_isb = x[IC(2, opt)];
				ppp_residual[idx].cmp_isb = x[IC(3, opt)];
				ppp_residual[idx].qzs_isb = x[IC(4, opt)];
				ppp_residual[idx].dion = dion_tec + dion;
				ppp_residual[idx].dtrp = dtrp;
				ppp_residual[idx].varerr = verr;

				
				ppp_residual[idx].residual = resid;
				
				ppp_residual[idx].rs[0] = (rs + i * 6)[0];
				ppp_residual[idx].rs[1] = (rs + i * 6)[1];
				ppp_residual[idx].rs[2] = (rs + i * 6)[2];
				res_cnt++;
			}

		}
		//constraint to external ionosphere correction
		/*if (opt->ionoopt == IONOOPT_CONST) {
			int nx_nv = nx * nv;
			int jj = II(sat, opt);
			v[nv] = dion_tec - x[jj];
			for (kk = 0; kk < nx; kk++) H[kk + nx_nv] = kk == jj ? 1.0 : 0.0;
			var[nv++] = 0.2;
		}*/

	}
	/* reject satellite with large and max post-fit residual */
	if (post && ne > 0) {
		vmax = ve[0]; maxobs = obsi[0]; maxfrq = frqi[0]; rej = 0;
		for (j = 1; j < ne; j++) {
			if (fabs(vmax) >= fabs(ve[j])) continue;
			vmax = ve[j]; maxobs = obsi[j]; maxfrq = frqi[j]; rej = j;
		}
		sat = obs[maxobs].sat;
		trace(2, "outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
			post, str, sat, maxfrq % 2 ? "P" : "L", maxfrq / 2 + 1, vmax, azel[1 + maxobs * 2] * R2D);
		exc[maxobs] = 1; rtk->ssat[sat - 1].rejc[maxfrq % 2]++; stat = 0;
		ve[rej] = 0;
	}
	for (i = 0; i < nv; i++) for (j = 0; j < nv; j++) {
		R[i + j * nv] = i == j ? var[i] : 0.0;
	}
	return post ? stat : nv;
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t* opt)
{
	return NX(opt);
}
/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t* rtk, const obsd_t* obs, int n, int stat)
{
	const prcopt_t* opt = &rtk->opt;
	int i, j;

	/* test # of valid satellites */
	rtk->sol.ns = 0;
	for (i = 0; i < n && i < MAXOBS; i++) {
		for (j = 0; j < opt->nf; j++) {
			if (!rtk->ssat[obs[i].sat - 1].vsat[j]) continue;
			rtk->ssat[obs[i].sat - 1].lock[j]++;
			rtk->ssat[obs[i].sat - 1].outc[j] = 0;
			if (j == 0) rtk->sol.ns++;
		}
	}
	int min_sat = MIN_NSAT_SOL;
	rtk->sol.stat = rtk->sol.ns < MIN_NSAT_SOL ? SOLQ_NONE : stat;

	if (rtk->sol.stat == SOLQ_FIX) {
		for (i = 0; i < 3; i++) {
			rtk->sol.rr[i] = rtk->xa[i];
			rtk->sol.qr[i] = (float)rtk->Pa[i + i * rtk->na];
		}
		rtk->sol.qr[3] = (float)rtk->Pa[1];
		rtk->sol.qr[4] = (float)rtk->Pa[1 + 2 * rtk->na];
		rtk->sol.qr[5] = (float)rtk->Pa[2];
	}
	else {
		for (i = 0; i < 3; i++) {
			rtk->sol.rr[i] = rtk->x[i];
			rtk->sol.qr[i] = (float)rtk->P[i + i * rtk->nx];
		}
		rtk->sol.qr[3] = (float)rtk->P[1];
		rtk->sol.qr[4] = (float)rtk->P[2 + rtk->nx];
		rtk->sol.qr[5] = (float)rtk->P[2];
	}
	
	rtk->sol.dtr[0] = rtk->x[IC(0, opt)];
	rtk->sol.dtr[1] = rtk->x[IC(1, opt)] - rtk->x[IC(0, opt)];
	rtk->sol.time = obs[0].time;

	for (i = 0; i < n && i < MAXOBS; i++) for (j = 0; j < opt->nf; j++) {
		rtk->ssat[obs[i].sat - 1].snr[j] = obs[i].SNR[j];
	}
	for (i = 0; i < MAXSAT; i++) for (j = 0; j < opt->nf; j++) {
		if (rtk->ssat[i].slip[j] & 3) rtk->ssat[i].slipc[j]++;
		if (rtk->ssat[i].fix[j] == 2 && stat != SOLQ_FIX) rtk->ssat[i].fix[j] = 1;
	}
}
/* test hold ambiguity -------------------------------------------------------*/
static int test_hold_amb(rtk_t* rtk)
{
	int i, j, stat = 0;

	/* no fix-and-hold mode */
	if (rtk->opt.modear != ARMODE_FIXHOLD) return 0;

	/* reset # of continuous fixed if new ambiguity introduced */
	for (i = 0; i < MAXSAT; i++) {
		if (rtk->ssat[i].fix[0] != 2 && rtk->ssat[i].fix[1] != 2) continue;
		for (j = 0; j < MAXSAT; j++) {
			if (rtk->ssat[j].fix[0] != 2 && rtk->ssat[j].fix[1] != 2) continue;
			if (!rtk->ambc[j].flags[i] || !rtk->ambc[i].flags[j]) stat = 1;
			rtk->ambc[j].flags[i] = rtk->ambc[i].flags[j] = 1;
		}
	}
	if (stat) {
		rtk->nfix = 0;
		return 0;
	}
	/* test # of continuous fixed */
	return ++rtk->nfix >= rtk->opt.minfix;
}

void print_ppp_residual(const residual_t* res, const int cnt) {
	for (int i = 0; i < cnt; i++) {
		fprintf(ppp_out_file, "%d %.5f %d %d %d %.2f %.2f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
			res[i].week, res[i].gpst, res[i].sys, res[i].prn, res[i].freq,
			res[i].el, res[i].snr,
			res[i].dtr, res[i].glo_isb, res[i].gal_isb,
			res[i].dion, res[i].dtrp,
			res[i].residual, res[i].varerr);
	}
}

/* precise point positioning -------------------------------------------------*/
extern void pppos(rtk_t* rtk, const obsd_t* obs, int n, const nav_t* nav)
{
	const prcopt_t* opt = &rtk->opt;
	double* rs, * dts, * var, * v, * H, * R, * azel, * xp, * Pp, dr[3] = { 0 }, std[3];
	char str[32];
	int i, j, nv, info, svh[MAXOBS], exc[MAXOBS] = { 0 }, stat = SOLQ_SINGLE;

	time2str(obs[0].time, str, 2);
	trace(3, "pppos   : time=%s nx=%d n=%d\n", str, rtk->nx, n);

	rs = mat(6, n); dts = mat(2, n); var = mat(1, n); azel = zeros(2, n);

	for (i = 0; i < MAXSAT; i++) for (j = 0; j < opt->nf; j++) rtk->ssat[i].fix[j] = 0;

	int satn = satno(0x04, 5);

	/* temporal update of ekf states */
	udstate_ppp(rtk, obs, n, nav);

	/* satellite positions and clocks */
	satposs(obs[0].time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);

	/* exclude measurements of eclipsing satellite (block IIA) */
	if (rtk->opt.posopt[3]) {
		testeclipse(obs, n, nav, rs);
	}
	/* earth tides correction */
	if (opt->tidecorr) {
		tidedisp(gpst2utc(obs[0].time), rtk->x, opt->tidecorr == 1 ? 1 : 7, &nav->erp,
			opt->odisp[0], dr);
	}
	nv = n * rtk->opt.nf * 2 + MAXSAT + 3;
	xp = mat(rtk->nx, 1); Pp = zeros(rtk->nx, rtk->nx);
	v = mat(nv, 1); H = mat(rtk->nx, nv); R = mat(nv, nv);

	for (i = 0; i < MAX_ITER; i++) {

		matcpy(xp, rtk->x, rtk->nx, 1);
		matcpy(Pp, rtk->P, rtk->nx, rtk->nx);

		/* prefit residuals */
		if (!(nv = ppp_res(0, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel))) {
			trace(2, "%s ppp (%d) no valid obs data\n", str, i + 1);
			break;
		}
		/* measurement update of ekf states */
		if ((info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {
			trace(2, "%s ppp (%d) filter error info=%d\n", str, i + 1, info);
			break;
		}
		/* postfit residuals */
		if (ppp_res(i + 1, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel)) {
			matcpy(rtk->x, xp, rtk->nx, 1);
			matcpy(rtk->P, Pp, rtk->nx, rtk->nx);
			stat = SOLQ_PPP;
			break;
		}
	}
	if (i >= MAX_ITER) {
		trace(2, "%s ppp (%d) iteration overflows\n", str, i);
	}
	if (stat == SOLQ_PPP) {

		/* ambiguity resolution in ppp */
		if (ppp_ar(rtk, obs, n, exc, nav, azel, xp, Pp) &&
			ppp_res(9, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel)) {

			matcpy(rtk->xa, xp, rtk->nx, 1);
			matcpy(rtk->Pa, Pp, rtk->nx, rtk->nx);
			
			for (i = 0; i < 3; i++) std[i] = sqrt(Pp[i + i * rtk->nx]);
			if (norm(std, 3) < MAX_STD_FIX) stat = SOLQ_FIX;
		}
		else {
			rtk->nfix = 0;
		}
		/* update solution status */
		update_stat(rtk, obs, n, stat);

		/* print residuals */
		print_ppp_residual(ppp_residual, nv);

		/* hold fixed ambiguities */
		if (stat == SOLQ_FIX && test_hold_amb(rtk)) {
			matcpy(rtk->x, xp, rtk->nx, 1);
			matcpy(rtk->P, Pp, rtk->nx, rtk->nx);
			trace(2, "%s hold ambiguity\n", str);
			rtk->nfix = 0;
		}
	}
	free(rs); free(dts); free(var); free(azel);
	free(xp); free(Pp); free(v); free(H); free(R);

}
