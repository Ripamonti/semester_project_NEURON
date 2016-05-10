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
#include<vector>
#include<problem.h>

class CN
{
  public:
    // Vectors to store the main elements of the linear system to be solved
    //   N.B. Even if we don't have a tridiagonal system the elements are stored
    //   as if we have it. This is for improving performance but we have to take it
    //   into account for the linear solver
    std::vector<double> upper_diagonal;
    std::vector<double> lower_diagonal;
    std::vector<double> mid_diagonal;
    std::vector<double> known_term;
    // Biological constants 
    static constexpr double c_m=1;
    static constexpr double g_l=0.0003;
    static constexpr double e_l=-54.3;
    static constexpr double R=180;
    static constexpr double a=0.5;
    static constexpr double g_K=0.036;
    static constexpr double e_K=-77;
    static constexpr double g_Na=0.12;
    static constexpr double e_Na=50;
    static constexpr double t_off=1;
    static constexpr double tau=1.74;
    static constexpr double g_syn=0.001;
    static constexpr double e_syn=0;
    // Store the total lengths of the three branches
    //   N.B. In this case we count for each case the branching node
    std::size_t n_L1;
    std::size_t n_L2;
    std::size_t n_L3;
    // Pointers to vector related to the gate states and the potential of the previous
    // stage
    std::vector<double>* n;	
    std::vector<double>* m;	
    std::vector<double>* h;
    std::vector<double>* y;
  public:
    // Variable used to have some sort of compatibility with the explicit method (ROCK2,RKC)
    // TODO: implement them in a useful way
    bool verbose;
    bool one_step;
    // Store the spatial step
    double s_step;
  public:
    // Constructor
    CN(CABLE* ode, bool onestep=true, bool verb=true);
    // Compute a step
    // TODO: implement some sort of control that can be checked through idid
    void advance(CABLE* ode, double& h, int& idid);
    // Build the main diagonal and the known term, that are the only elements of the
    // matrix that change through the simulation
    void build_diag(CABLE* ode, double dt);
    // Modified linear solver for quasi-tridiagonal system
    void modified_Thomas();
};

