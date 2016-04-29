#ifndef PROBLEM_H
#define	PROBLEM_H

#include <vector>
#include "Ode.h"
//#include "Ode.cpp"

class CABLE : public ODE
{
  public:
    // Pointers to the gate state values
    std::vector<double>* n;	
    std::vector<double>* m;	
    std::vector<double>* h;	
    int n_L1,n_L2,n_L3;
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
    
  public:
    CABLE(Mesh& mesh,bool intrho,double t,std::vector<double>& init);	// Constructor
    void get_gate_state(std::vector<double>& sodium,std::vector<double>& potassium,std::vector<double>& leak);	// Get a pointer to the vectors containing the gates states
    void set_Dirichlet(double t, std::vector<double>& y, bool use=false);
    void rhs(double t, std::vector<double>& y, std::vector<double>& f);
};

class GATE_N : public ODE
{
  public:
    std::vector<double>* v;	// Pointer to the potential value
    double alpha_n,beta_n;	// just to make the code clear
  public:
    GATE_N(Mesh& mesh,bool intrho,double t,std::vector<double>& init);
    void get_potential(std::vector<double>& y);	// Get a pointer to the vector containing the potential
    void set_Dirichlet(double t, std::vector<double>& y, bool use=false);
    void rhs(double t, std::vector<double>& y, std::vector<double>& f);
};

class GATE_M : public ODE
{
  public:
    std::vector<double>* v;	// Pointer to the potential value
    double alpha_m,beta_m;	// just to make the code clear
  public:
    GATE_M(Mesh& mesh,bool intrho,double t,std::vector<double>& init);
    void get_potential(std::vector<double>& y);	// Get a pointer to the vector containing the potential
    void set_Dirichlet(double t, std::vector<double>& y, bool use=false);
    void rhs(double t, std::vector<double>& y, std::vector<double>& f);
};

class GATE_H : public ODE
{
  public:
    std::vector<double>* v;	// Pointer to the potential value
    double alpha_h,beta_h;	// just to make the code clear
  public:
    GATE_H(Mesh& mesh,bool intrho,double t,std::vector<double>& init);
    void get_potential(std::vector<double>& y);	// Get a pointer to the vector containing the potential
    void set_Dirichlet(double t, std::vector<double>& y, bool use=false);
    void rhs(double t, std::vector<double>& y, std::vector<double>& f);
};

#endif	/* PROBLEM_H */
