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

#ifndef ODE_H
#define	ODE_H

#include "Mesh.h"
#include <vector>

class ODE
{
  public:
    std::vector<double> un;		// Initial value

    // Allocate vector that keep in memory the step modified by RK
    std::vector<double> yn;		// Space for initial value at each step
    std::vector<double> ynpu;		// Space for the next computed solution in the step

    double time;			// Actual time in integration
    double tend;			// Ending time of computation

    bool cte_rho;			// TRUE if the spectrum is constant, FALSE otherwise
    double s_step;			// Spatial step size

    // Integers used in the computation of the rhs to store the positions of the boundaries and the branching node
    int limit_branch_A;
    int limit_branch_B;
    int limit_branch_C;
    
    std::vector<double> fn;

    // Support vectors often used in the code to keep track of intermediate steps. They are declared here and initialized in the defined in the compiled file
    // to save some time during computation
    std::vector<double> dz1;
    std::vector<double> z;
    std::vector<double> fnpu;
    std::vector<double> tmp1;
 
    ODE(Mesh& mesh,bool intrho,double t);
    void print_info();  
    std::vector<double>& get_un();  
    void set_un(std::vector<double>&); 
    void rho(double& eigmax); 
    void rhs(double t, std::vector<double>& y, std::vector<double>& f);
    void set_Dirichlet(double t, std::vector<double>& y, bool use=false);
    double normalized_L2_norm(std::vector<double>& v);
   // ~ODE();
};

#endif	/* MESH_H */
