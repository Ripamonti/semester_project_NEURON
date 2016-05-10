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
    // Initial value
    //   N.B. During the first step this stores the initial value.
    //   Later in the simulation it will be used to store the value of the previous step
    std::vector<double> un;
    // Allocate vector that keep in memory the step modified by RK
    std::vector<double> yn;		// Space for initial value at each step
    std::vector<double> ynpu;		// Space for the next computed solution in the step
    // Time variables
    double time;			// Actual time in integration
    double tend;			// Ending time of computation
 
    bool cte_rho;			// TRUE if the spectrum is constant, FALSE otherwise
    double s_step;			// Spatial step size

    // Integers used in the computation of the rhs to store the positions of the boundaries and the branching node
    int limit_branch_A;
    int limit_branch_B;
    int limit_branch_C;
    // Support vector used during the computation of the spectral radius through the power method
    std::vector<double> fn;
    // Support vectors often used in the code to keep track of intermediate steps. They are declared here
    std::vector<double> dz1;
    std::vector<double> z;
    std::vector<double> fnpu;
    std::vector<double> tmp1;
 public:
    // Constructor
    ODE(Mesh& mesh,bool intrho,double t,std::vector<double>& init);
    // Print info regarding the problem
    void print_info();  
    // Method to get the initial value/ previous step solution
    //  N.B. A method is defined for compatibility with the ROCK2/RKC solver
    std::vector<double>& get_un();  
    // Method to set the vector un
    //  N.B. A method is defined for compatibility with the ROCK2/RKC solver
    void set_un(std::vector<double>&); 
    // Method used to compute internally the spectral radius
    void rho(double& eigmax); 
    // Method to define the rhs
    virtual void rhs(double t, std::vector<double>& y, std::vector<double>& f)=0;
    // Method used to impose Dirichlet boundary condition if needed
    virtual void set_Dirichlet(double t, std::vector<double>& y, bool use=false)=0;
    //  N.B. The previous two methods are implemented as virtual functions because
    //       are problem dependent: for different problems one has to define manually
    //       their behaviour in the derived class
    // Compute the normalize L2 norm
    double normalized_L2_norm(std::vector<double>& v);
    // TODO: implement the destructor of the class (I don't know if the standard is enough)
    // ~ODE();
};

#endif	/* MESH_H */
