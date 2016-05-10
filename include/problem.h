#ifndef PROBLEM_H
#define	PROBLEM_H

#include <vector>
#include "Ode.h"

// Definition and declaration of derived classes for potential equation
class CABLE : public ODE
{
  public:
    // Pointers to the gate state values
    std::vector<double>* n;	
    std::vector<double>* m;	
    std::vector<double>* h;	
    // Variables used to store the length of each branch
    std::size_t n_L1;
    std::size_t n_L2;
    std::size_t n_L3;
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
  public:
    // Constructor
    CABLE(Mesh& mesh,bool intrho,double t,std::vector<double>& init);
    // Method to get the states (used in rhs computation)
    void get_gate_state(std::vector<double>& sodium,std::vector<double>& potassium,std::vector<double>& leak);
    // Method to set the Dirichlet conditions 
    void set_Dirichlet(double t, std::vector<double>& y, bool use=false);
    // Method to define the right hand side
    void rhs(double t, std::vector<double>& y, std::vector<double>& f);
};

// Definition and declaration of derived class for gate n equation
class GATE_N : public ODE
{
  public:
    // Pointer to the potential value
    std::vector<double>* v;
    // TODO: delete these variables... They are not used now!!	
    double alpha_n,beta_n;
  public:
    // Constructor
    GATE_N(Mesh& mesh,bool intrho,double t,std::vector<double>& init);
    // Method to get the potential (used in rhs computation)
    void get_potential(std::vector<double>& y);
    // Method to set the Dirichlet conditions  
    void set_Dirichlet(double t, std::vector<double>& y, bool use=false);
    // Method to define the right hand side
    void rhs(double t, std::vector<double>& y, std::vector<double>& f);
};

// Definition and declaration of derived class for gate m equation
class GATE_M : public ODE
{
  public:
    // Pointer to the potential value
    std::vector<double>* v;
    // TODO: delete these variables... They are not used now!!
    double alpha_m,beta_m;
  public:
    // Constructor
    GATE_M(Mesh& mesh,bool intrho,double t,std::vector<double>& init);
    // Method to get the potential (used in rhs computation)
    void get_potential(std::vector<double>& y);	
    // Method to set the Dirichlet conditions  
    void set_Dirichlet(double t, std::vector<double>& y, bool use=false);
    // Method to define the right hand side
    void rhs(double t, std::vector<double>& y, std::vector<double>& f);
};

// Definition and declaration of derived class for gate h equation
class GATE_H : public ODE
{
  public:
    // Pointer to the potential value
    std::vector<double>* v;	
    // TODO: delete these variables... They are not used now!!
    double alpha_h,beta_h;	
  public:
    // Constructor
    GATE_H(Mesh& mesh,bool intrho,double t,std::vector<double>& init);
    // Method to get the potential (used in rhs computation)
    void get_potential(std::vector<double>& y);
    // Method to set the Dirichlet conditions  
    void set_Dirichlet(double t, std::vector<double>& y, bool use=false);
    // Method to define the right hand side
    void rhs(double t, std::vector<double>& y, std::vector<double>& f);
};

#endif	/* PROBLEM_H */
