#include "rtklib.h"

#define SQR(x)      ((x)*(x))

void buildDatPos(gtime_t t, const export_t* data);

/* get group delay parameter (m) ---------------------------------------------*/
static double gettaugd(int sat, const nav_t* nav, int type)
{
    int i, sys = satsys(sat, NULL);

    if (sys == SYS_GLO) {
        for (i = 0; i < nav->ng; i++) {
            if (nav->geph[i].sat == sat) break;
        }
        return (i >= nav->ng) ? 0.0 : -nav->geph[i].dtaun * CLIGHT;
    }
    else {
        for (i = 0; i < nav->n; i++) {
            if (nav->eph[i].sat == sat) break;
        }
        return (i >= nav->n) ? 0.0 : nav->eph[i].tgd[type] * CLIGHT;
    }
}

void allocate(export_t* data, int n) {

    // Allocate dynamic memmory
    data->sat = mat(1, n);
    data->geom = mat(1, n);
    data->snr = mat(2, n);
    data->range = mat(2, n);
    data->phase = mat(2, n);
    data->doppler = mat(2, n);

}

void initdata(export_t* data, int n) {

    data->ntot = 0;
    data->azel[0] = 0.0;
    data->azel[1] = 0.0;
    data->rcvclk = 0.0;

    memset(data->range, 0, sizeof(double) * 2 * n);
    memset(data->phase, 0, sizeof(double) * 2 * n);
    memset(data->doppler, 0, sizeof(double) * 2 * n);
    memset(data->snr, 0, sizeof(double) * 2 * n);
    memset(data->sat, 0, sizeof(double) * n);
    memset(data->geom, 0, sizeof(double) * n);
}


void preparedata(export_t* data, const obsd_t* obs, const nav_t* nav, int i, int* ni, 
    int sat, int sys, double freq, double dionL1Gps, double dtrp, const double* azel, double r) {

    int ind = (*ni);

    // Prepare data for export
    double taugd_corr1 = 0.0;
    double taugd_corr2 = 0.0;

    // Extract taugd for L1 and L5 frequency
    if (sys == SYS_GPS || sys == SYS_QZS) { /* L1 */
        taugd_corr1 = gettaugd(sat, nav, 0); /* TGD (m) */
        if ((freq = sat2freq(sat, obs[i].code[2], nav)))
            taugd_corr2 = taugd_corr1 * SQR(FREQ1 / freq); /* L5 */
    }
    else if (sys == SYS_GLO) { /* G1 */
        double gamma = SQR(FREQ1_GLO / FREQ2_GLO);
        double b1 = gettaugd(sat, nav, 0); /* -dtaun (m) */
        taugd_corr1 = b1 / (gamma - 1.0);
    }
    else if (sys == SYS_GAL) { /* E1 */
        if (getseleph(SYS_GAL)) taugd_corr1 = gettaugd(sat, nav, 0); /* BGD_E1E5a */
        else                    taugd_corr1 = gettaugd(sat, nav, 1); /* BGD_E1E5b */

        if ((freq = sat2freq(sat, obs[i].code[2], nav)))
            taugd_corr2 = taugd_corr1 * SQR(FREQ1 / freq); /* E5 */

    }
    else if (sys == SYS_CMP) { /* B1I/B1Cp/B1Cd */
        if (obs->code[0] == CODE_L2I) taugd_corr1 = gettaugd(sat, nav, 0); /* TGD_B1I */
        else if (obs->code[0] == CODE_L1P) taugd_corr1 = gettaugd(sat, nav, 2); /* TGD_B1Cp */
        else taugd_corr1 = gettaugd(sat, nav, 2) + gettaugd(sat, nav, 4); /* TGD_B1Cp+ISC_B1Cd */
    }

    double dion1 = 0.0;
    double dion2 = 0.0;

    // Scale Iono correction
    if ((freq = sat2freq(sat, obs[i].code[0], nav)))
        dion1 = dionL1Gps * SQR(FREQ1 / freq);
    if ((freq = sat2freq(sat, obs[i].code[2], nav)))
        dion2 = dionL1Gps * SQR(FREQ1 / freq);

    double* p_range = data->range + i*2;
    if (obs[i].P[0] != 0.0)
        p_range[0] = obs[i].P[0] - dion1 - dtrp - taugd_corr1;
    else
        p_range[0] = 0.0;
    
    if (obs[i].P[2] != 0.0)
        p_range[1] = obs[i].P[2] - dion2 - dtrp - taugd_corr2;
    else
        p_range[1] = 0.0;

    double* p_snr = data->snr + i * 2;
    if (obs[i].SNR[0] * SNR_UNIT != 0.0)
        p_snr[0] = obs[i].SNR[0] * SNR_UNIT;
    else
        p_snr[0] = 0.0;

    if (obs[i].SNR[2] * SNR_UNIT != 0.0)
        p_snr[1] = obs[i].SNR[1] * SNR_UNIT;
    else
        p_snr[1] = 0.0;

    double* p_phase = data->phase + i * 2;
    if (obs[i].L[0] != 0.0)
        p_phase[0] = obs[i].L[0] + dion1 - dtrp;
    else
        p_phase[0] = 0.0;

    if (obs[i].L[2] != 0.0)
        p_phase[1] = obs[i].L[2] + dion2 - dtrp;
    else
        p_phase[1] = 0.0;

    double* p_doppler = data->doppler + i*2;
    p_doppler[0] = obs[i].D[0];
    p_doppler[1] = obs[i].D[2];

    data->sat[ind] = sat;
    data->geom[ind] = r;

    data->azel[0] = (azel + i*2)[0];
    data->azel[1] = (azel + i*2)[1];

    (*ni)++;
}

void senddata(export_t* data, gtime_t time, const double* rs, const double* dts, double dtr, int ni) {

    data->rcvclk = dtr;
    data->ntot = ni;
    data->satpos = rs;
    data->satclk = dts;
    data->rcvclk = dtr;

   // buildDatPos(time, data);
}


void clearmemory(const export_t* data) {

    /* Free allocated memory */
    free(data->range);
    free(data->phase);
    free(data->doppler);
    free(data->sat);
    free(data->snr);
    free(data->geom);

}