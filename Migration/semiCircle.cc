/*
 * semiCircle.cc - 3D semi-circle migration for GPR data
 * by Morten Harms 2026
 *
 * Compile from the MATLAB command window:
 *   mex CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' semiCircle.cc
 *
 * Inputs:
 *   0  cube            double Nx x Ny x ns  - 3D array of trace data
 *   1  Besetzung       double Nx x Ny       - occupied bin map
 *   2  topoMig         double Nx x Ny       - surface topography (depth below reference)
 *   3  Nx              scalar int           - number of x bins
 *   4  Ny              scalar int           - number of y bins
 *   5  Nz              scalar int           - number of z bins in output
 *   6  ns              scalar int           - number of time samples per trace
 *   7  nTraces         scalar int           - total number of occupied traces
 *   8  dz              scalar double        - depth sampling interval [m]
 *   9  alpha0          scalar double        - maximum aperture angle [degrees]
 *  10  xaxis           double 1 x Nx        - x coordinates
 *  11  yaxis           double 1 x Ny        - y coordinates
 *  12  zaxis           double 1 x Nz        - depth axis
 *  13  v               scalar double        - migration velocity [m/ns]
 *  14  dt              scalar double        - time sampling interval [ns]
 *  15  maskIn          double Nx x Ny       - valid migration bin mask
 *  16  nCores          scalar int           - number of OpenMP threads requested
 *
 * Output:
 *   migcube  double Nx x Ny x Nz  - migrated 3D data cube
 */

#include "mex.h"
#include <omp.h>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 17)
        mexErrMsgIdAndTxt("semiCircle:nrhs", "17 input arguments required.");
    if (nlhs > 1)
        mexErrMsgIdAndTxt("semiCircle:nlhs", "At most 1 output argument.");

    /* --- set number of OpenMP threads ---
     * Clamp requested thread count to the number of available cores.  */
    int nCores     = static_cast<int>(mxGetScalar(prhs[16]));
    const int nMax = omp_get_num_procs();  // actual hardware processor count
    if (nCores > nMax) {
        nCores = nMax;
    }
    omp_set_num_threads(nCores);
    mexPrintf("OpenMP with %d of %d threads\n", nCores,nMax);

    /* --- scalar inputs --- */
    const int    Nx       = static_cast<int>(mxGetScalar(prhs[3]));
    const int    Ny       = static_cast<int>(mxGetScalar(prhs[4]));
    const int    Nz       = static_cast<int>(mxGetScalar(prhs[5]));
    const int    ns       = static_cast<int>(mxGetScalar(prhs[6]));
    const int    nTraces  = static_cast<int>(mxGetScalar(prhs[7]));
    const double dz       = mxGetScalar(prhs[8]);
    const double alpha0   = mxGetScalar(prhs[9]);
    const double v        = mxGetScalar(prhs[13]);
    const double dt       = mxGetScalar(prhs[14]);

    /* --- pointers to input arrays --- */
    const double *mxDatacube  = mxGetPr(prhs[0]);
    const double *mxBesetzung = mxGetPr(prhs[1]);
    const double *mxTopo      = mxGetPr(prhs[2]);
    const double *mxXaxis     = mxGetPr(prhs[10]);
    const double *mxYaxis     = mxGetPr(prhs[11]);
    const double *mxZaxis     = mxGetPr(prhs[12]);
    const double *mxMaskIn    = mxGetPr(prhs[15]);

    /* --- dimension shortcuts --- */
    const int NxNy   = Nx * Ny;
    const int NxNyNz = Nx * Ny * Nz;

    /* --- copy MATLAB column-major inputs into flat C row-major working arrays ---
     * Row-major layout keeps the innermost loop index (depth/time) contiguous
     * in memory, which is cache-efficient for the migration loop below.       */
    std::vector<double> datacube  (Nx * Ny * ns, 0.0);
    std::vector<double> besetzung (NxNy, 0.0);
    std::vector<double> topography(NxNy, 0.0);
    std::vector<double> mask_in   (NxNy, 0.0);
    std::vector<double> x_axis    (Nx,   0.0);
    std::vector<double> y_axis    (Ny,   0.0);
    std::vector<double> z_axis    (Nz,   0.0);

    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++) {
            besetzung [i*Ny + j] = mxBesetzung[i + Nx*j];
            topography[i*Ny + j] = mxTopo     [i + Nx*j];
            mask_in   [i*Ny + j] = mxMaskIn   [i + Nx*j];
        }

    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < ns; k++)
                datacube[i*Ny*ns + j*ns + k] = mxDatacube[i + Nx*j + NxNy*k];

    for (int i = 0; i < Nx; i++) x_axis[i] = mxXaxis[i];
    for (int j = 0; j < Ny; j++) y_axis[j] = mxYaxis[j];
    for (int k = 0; k < Nz; k++) z_axis[k] = mxZaxis[k];

    /* --- output and accumulation arrays --- */
    std::vector<double> migcube  (NxNyNz, 0.0);
    std::vector<double> countcube(NxNyNz, 0.0);
    std::vector<double> mask_buf (NxNyNz, 0.0);

    /* --- pre-migration constants --- */
    const double alpha0Rad    = (alpha0 / 360.0) * (2.0 * M_PI);
    const double halfAperture = alpha0Rad / 2.0;           // aperture half-angle [rad]
    const double b            = (2.0 * M_PI) / alpha0Rad; // taper frequency factor

    /* --- print parameters and start time --- */
    char text[100];
    std::time_t now = std::time(nullptr);
    struct std::tm *t = std::localtime(&now);
    mexPrintf("alpha0 = %f v = %f dt = %f dz = %f\n", alpha0, v, dt, dz);
    std::strftime(text, sizeof(text)-1, "%d %m %Y %H:%M", t);
    mexPrintf("Start date: %s\n", text);

    int tr = 0;
    mexPrintf("- migrating trace 0 of %d (0%%)\n", nTraces);
    mexEvalString("drawnow");

    /* --- main migration loop ---
     *
     * For each occupied antenna position (ixt, iyt):
     *   For each valid migration bin (ixd, iyd, iz):
     *     - Compute the ray path length from antenna to bin
     *     - Look up the corresponding time sample in the input trace
     *     - Apply an angle-dependent cosine taper
     *     - Accumulate into the migration output cube               */
    for (int ixt = 0; ixt < Nx; ixt++)
    {
        for (int iyt = 0; iyt < Ny; iyt++)
        {
            if (besetzung[ixt*Ny + iyt] > 0)
            {
                tr++;
                if (tr % 100 == 0 || tr == nTraces) {
                    mexPrintf("- migrating trace %d of %d (%.1f%%)\n",
                              tr, nTraces, 100.0 * tr / nTraces);
                    mexEvalString("drawnow");
                }

                // antenna position at the surface
                const double xAnt = x_axis[ixt];
                const double yAnt = y_axis[iyt];
                const double zAnt = topography[ixt*Ny + iyt];

                std::fill(mask_buf.begin(), mask_buf.end(), 0.0);

                #pragma omp parallel for collapse(2)
                for (int ixd = 0; ixd < Nx; ixd++)
                {
                    for (int iyd = 0; iyd < Ny; iyd++)
                    {
                        if (mask_in[ixd*Ny + iyd] > 0)
                        {
                            // vector (d): from antenna to migration bin projected onto surface
                            const double xMig  = x_axis[ixd] - xAnt;
                            const double xMig2 = xMig * xMig;
                            const double yMig  = y_axis[iyd] - yAnt;
                            const double yMig2 = yMig * yMig;
                            const double zSurf = topography[ixd*Ny + iyd] - zAnt;
                            const double lengthD = std::sqrt(xMig2 + yMig2 + zSurf*zSurf);

                            const int iz_start = static_cast<int>(
                                std::round(topography[ixd*Ny + iyd] / dz));

                            for (int iz = iz_start; iz < Nz; iz++)
                            {
                                // vector (r): ray path from antenna to migration bin
                                const double zMig    = z_axis[iz] - zAnt;
                                const double lengthR = std::sqrt(xMig2 + yMig2 + zMig*zMig);

                                // angle between r and d, measured from vertical [rad]
                                // (0 directly below antenna, pi/2 at surface level)
                                double alphaRad;
                                if (ixd == ixt && iyd == iyt) {
                                    alphaRad = 0.0;
                                } else {
                                    // dot product of r and d: r·d = xMig² + yMig² + zMig*zSurf
                                    const double skalarprod = xMig2 + yMig2 + zMig * zSurf;
                                    const double cosVal     = skalarprod / (lengthR * lengthD);
                                    // clamp to [-1,1] to guard acos against floating-point rounding
                                    const double cosC = std::max(-1.0, std::min(1.0, cosVal));
                                    alphaRad = (M_PI / 2.0) - std::acos(cosC);
                                }

                                // two-way travel time to input sample index
                                const double twt  = lengthR * 2.0 / v;
                                const int    indx = static_cast<int>(std::round(twt / dt));

                                if (indx >= 0 && indx < ns && alphaRad < halfAperture)
                                {
                                    // cosine taper: 1 at beam centre, 0 at aperture edge
                                    const double taper  = 1.0 - ((std::cos(alphaRad * b) + 1.0) * 0.5);
                                    const double sample = datacube[ixt*Ny*ns + iyt*ns + indx];
                                    const int    midx   = ixd*Ny*Nz + iyd*Nz + iz;

                                    mask_buf[midx]  += sample * taper;
                                    countcube[midx] += 1.0;
                                }
                            } // iz
                        } // mask_in
                    } // iyd
                } // ixd

                for (int idx = 0; idx < NxNyNz; idx++)
                    migcube[idx] += mask_buf[idx];

            } // besetzung
        } // iyt
    } // ixt

    mexPrintf("finished! (%d traces migrated)\n\n", tr);
    mexEvalString("drawnow");

    /* --- normalise each bin by its total hit count --- */
    for (int idx = 0; idx < NxNyNz; idx++)
        if (countcube[idx] > 0.0)
            migcube[idx] /= countcube[idx];

    /* --- copy result into MATLAB output array ---
     * Internal arrays are row-major; MATLAB expects column-major,
     * so indices are transposed on the way out.                    */
    mwSize dims[3] = { (mwSize)Nx, (mwSize)Ny, (mwSize)Nz };
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *result = mxGetPr(plhs[0]);

    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
                result[i + Nx*j + NxNy*k] = migcube[i*Ny*Nz + j*Nz + k];

    now = std::time(nullptr);
    t   = std::localtime(&now);
    std::strftime(text, sizeof(text)-1, "%d %m %Y %H:%M", t);
    mexPrintf("alpha0 = %f v = %f dt = %f dz = %f  %s\n", alpha0, v, dt, dz, text);
}
