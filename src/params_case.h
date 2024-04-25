/*
 * File Name: params_case.h
 * Description: Declare CaseParams class used set parameters by users
 * Created by Hao-Xin on 2021/10/15.
 *
 */


#ifndef HUBBARD_SRC_PARAMS_CASE_H_
#define HUBBARD_SRC_PARAMS_CASE_H_

#include "qlmps/case_params_parser.h"
using qlmps::CaseParamsParserBasic;

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Lx = ParseInt("Lx");
    Numhole = ParseInt("Numhole");
    t = ParseDouble("t");
    t_perp = ParseDouble("t_perp");
    t2 = ParseDouble("t2");
    U = ParseDouble("U");
    V = ParseDouble("V");
    V2 = ParseDouble("V2");
    J = ParseDouble("J");
    Sweeps = ParseInt("Sweeps");
    tau = ParseDouble("tau");
    steps = ParseInt("steps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
    Threads = ParseInt("Threads");
    Perturbation=ParseDouble("Perturbation");
    wavelength = ParseInt("wavelength");
    noise = ParseDoubleVec("noise");
    CorrelationMode = ParseInt("CorrelationMode");
  }

  size_t Lx;           //System size
  size_t Numhole;     //number of holes
  double t;           //hopping
  double t_perp;
  double t2;          //NNN hopping
  double U;           //on-site potential
  double V;           //nearest neighbor interaction
  double V2;          //NNN interaction
  double J;           //hopping square terms
  size_t Sweeps;      //sweep times for DMRG
  double tau;         //tau for TDVP
  size_t steps;       //steps for TDVP, note tau * steps = total time
  size_t Dmin;        //bond dimension
  size_t Dmax;        //bond dimension
  double CutOff;      //for SVD truncation
  double LanczErr;    //lanczos error
  size_t MaxLanczIter;//maximum lanczos iteration times
  size_t Threads;     //the number of thread
  double Perturbation;//charge density perturbation amplitude
  size_t wavelength;  //wave length of charge density perturbation
  std::vector<double> noise;
  size_t CorrelationMode; //0 : spin up, <c^dag(t) c(0)>
                          //1 : spin up, <c(t) c^dag(0)>
                          //2 : spin down, <c^dag(t) c(0)>
                          //3 : spin down, <c(t) c^dag(0)>
};

#endif //HUBBARD_SRC_PARAMS_CASE_H_
