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

#ifndef OPERATION_H
#define OPERATION_H

#include <algorithm>
#include <numeric>

// Function to sum two vector: v=v+c*w
void add(std::vector<double>& v, double c, std::vector<double>& w);
// function to add a constant to the vector: v=v+c*(unitary vector)
void add(std::vector<double>& v, double c);
// Function to sum two vectors: v=v+w (More efficient for this simpler case)
void add(std::vector<double>& v, std::vector<double>& w);
// Multiply a vector for a constatn: v=c*v
void scale(std::vector<double>& v, double c);
// Compute the L2 norm
double l2_norm(std::vector<double> const& u);
// Compute the discrete L2 norm
double l2_cont_norm(std::vector<double> const& u, double c);

#endif	/* OPERATION_H */
