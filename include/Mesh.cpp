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

#include "Mesh.h"
#include <cmath>

Mesh::Mesh(double dx): h_space(dx)
{
// Compute the number of the elements of the mesh given dx 
  n_elem = round(0.150/dx+2*0.250/dx+4);
  n_L1 = round(0.150/dx+2.);
  n_L2 = round(0.250/dx+2.);
  n_L3 = round(0.250/dx+2.);
// Reserve space in memory for the grid
  grid.resize(n_elem);
// Fill the points in the grid paying attention to the boundaries, where we have dx/2
// As center of the branch we the branching node
// First branch
  grid[0]=(n_L1-2)*h_space;
  grid[1]=grid[0]-dx/2;
  for (std::size_t i=2;i<(n_L1-1);i++)
  {
    grid[i]=grid[i-1]-dx;
  }  
// Second branch
  grid[n_L1-1]=(n_L2-2)*h_space;
  grid[n_L1]=grid[n_L1-1]-dx/2;
  for (std::size_t i=n_L1+1;i<(n_L1+n_L2-2);i++)
  {
    grid[i]=grid[i-1]-dx;
  }
// Third branch
  grid[n_L1+n_L2-2]=0;
  grid[n_L1+n_L2-1]=dx/2;
  for (std::size_t i=n_L1+n_L2; i<=(n_L1+n_L2+n_L3-4);i++)
  {
    grid[i]=grid[i-1]+dx;
  }
  grid[n_L1+n_L2+n_L3-3]=grid[n_L1+n_L2+n_L3-4]+dx/2;
}

void Mesh::print_info()
{
  printf("\n----------   Mesh informations   ------------------\n");
  printf("Number of elements: %zu\n",n_elem);
  printf("Grid space: %f\n",h_space);
  printf("----------------------------------------------------\n");
}

void Mesh::print_all()
{
  printf("List of all the nodes %lu \n",grid.size());
  for (std::size_t i=0; i<grid.size(); i++)
    printf(" %f ",grid[i]);
  fflush(stdout);
}

