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

#include "Ode.h"
#include <cmath>
#include <numeric>

ODE::ODE(Mesh& mesh,bool intrho,double t,std::vector<double>& init) : un(init), time(0) , tend(t), cte_rho(false), s_step(mesh.h_space), limit_branch_A(mesh.n_L1), limit_branch_B(mesh.n_L2)
, limit_branch_C(mesh.n_L3){
  // Reserve space for auxiliary variables used during the time step
  yn.resize(mesh.n_elem);
  ynpu.resize(mesh.n_elem);

  // Reserve space for auxiliary variables used during spectral radius computation
  dz1.resize(mesh.n_elem);
  z.resize(mesh.n_elem);
  fnpu.resize(mesh.n_elem);
  tmp1.resize(mesh.n_elem); // Probailmente da togliere..

  // Space for the evaluation
  fn.resize(mesh.n_elem);
}

void ODE::print_info()
{
  // TODO
}

std::vector<double>& ODE::get_un()
{
  return un;
}

void ODE::set_un(std::vector<double>& v)
{
  un=v;
}

void ODE::rho(double& eigmax)
{
  eigmax= 4/(s_step*s_step);
}

double ODE::normalized_L2_norm(std::vector<double>& u)
{
  return sqrt(s_step*std::inner_product(u.begin(), u.end(),u.begin(),0.0));
}
/*
ODE::~ODE(void)
{
  
}
*/
