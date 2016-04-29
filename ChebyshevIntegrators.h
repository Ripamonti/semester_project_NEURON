/*
 * © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
 * Laboratory ANMC, 2016
 * See the LICENSE.TXT file for more details.
 * 
 * Authors
 * The original Fortran codes were developed by Assyr Abdulle (for ROCK2) and 
 * Sommeijer, Shampine, Verwer (for RKC). Translation to C++ by Giacomo Rosilho
 * de Souza.
 * 
 * The ROCK2 code is described in 
 * 
 * Abdulle, Assyr, and Alexei A. Medovikov. 
 * "Second order Chebyshev methods based on orthogonal polynomials." 
 * Numerische Mathematik 90.1 (2001): 1-18.	
 * 
 * and the RKC code in
 * 
 * Sommeijer, B. P., L. F. Shampine, and J. G. Verwer. 
 * "RKC: An explicit solver for parabolic PDEs." 
 * Journal of Computational and Applied Mathematics 88.2 (1998): 315-326.
 * 
 * Please cite these articles in any publication describing research
 * performed using the software.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CHEBYSHEVINTEGRATORS_H
#define	CHEBYSHEVINTEGRATORS_H

#include "TimeIntegrator.h"
#include "Ode.h"
#include "operation.h"

class RKC: public TimeIntegrator<ODE>
{
public:
    RKC(bool onestep=true, bool verb=true, bool dtadap=true, 
        double atol=1e-2, double rtol=1e-2, bool intrho=false, bool scalartol=true);
protected:

    void rtstep(ODE* ode, const double t, const double& h, vector<double>& y,
                   vector<double>& yn);
    
    void update_n_stages(double& h);
    void init_coeffs(double *w, int s, double *bj, double *thj, double *zj, double *dzj, 
                     double *d2zj, double& kappa);
    void update_coeffs(double *w, double *bj, double *zj, double *dzj, double *d2zj, 
                       double& mu, double& nu, double& kappa, double& ajm1, double* thj);
    void shift_coeffs(double *bj, double *zj, double *dzj, double *d2zj, double* thj);
};

class ROCK2: public TimeIntegrator<ODE>
{
public:
    ROCK2(bool onestep=true, bool verb=true, bool dtadap=true, 
          double atol=1e-2, double rtol=1e-2, bool intrho=false, bool scalartol=true);
    
protected:

    void rtstep(ODE* ode, const double t, const double& h, vector<double>& y,
                   vector<double>& yn);

    void update_n_stages(double& h);
    void mdegr(int& mdeg, int mp[]);
    
 
protected:
    int mp[2];  ///<It is used in order to find the algorithm precomputed coefficients in the tables.
    
    static int ms[46];      ///<Array of coefficients.
    static double fp1[46];    ///<Array of coefficients.
    static double fp2[46];    ///<Array of coefficients.
    static double recf[4476]; ///<Array of coefficients.
};

#endif	/* CHEBYSHEVINTEGRATORS_H */

