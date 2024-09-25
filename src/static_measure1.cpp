/*
 * File Name: static_measure1.cpp
 * Description: measure static one point function
 * Created by Hao-Xin on 2021/11/09.
 *
 */

#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <ctime>
#include "qldouble.h"
#include "operators.h"
#include "params_case.h"
#include "myutil.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);
  size_t Lx = params.Lx;
  size_t N = 2 * Lx;

  clock_t startTime,endTime;
  startTime = clock();

  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);


  using namespace plain_hubbard;
  OperatorInitial();
  const SiteVec<TenElemT, U1U1QN> sites=SiteVec<TenElemT, U1U1QN>(N, pb_outF);


  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
  FiniteMPST mps(sites);
  mps.Load();
  cout << "mps loaded" <<endl;


  cout << "bond dimension of middle mps = " ;
  cout << mps[N/2].GetShape()[0] <<endl;
  Timer one_site_timer("measure one site operators");
  MeasureOneSiteOp(mps, {sz, nf}, {"static_sz","static_nf"});
  cout << "measured one point function.<====" <<endl;
  one_site_timer.PrintElapsed();

  endTime = clock();
  cout << "CPU Time : " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}