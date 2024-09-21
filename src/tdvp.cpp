/*
 * File Name: tdvp.cpp
 * Description: Calculate the single particle correlation for 1D Hubbard model
 * Created by Hao-Xin on 2021/11/07.
 *
 */

#include "qlten/qlten.h"
#include "qlmps/qlmps.h"
#include <time.h>
#include "params_case.h"
#include "qldouble.h"
#include "myutil.h"
#include "operators.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env(mpi::threading::multiple);
  mpi::communicator world;

  CaseParams params(argv[1]);
  size_t Lx = params.Lx;
  size_t N = 2 * Lx;

  double t = params.t, U = params.U, V = params.V;
  using std::cout;
  if (world.rank() == kMasterRank) {
    cout << " ****** 1D Hubbard Model Parameter List ****** " << "\n";
    cout << "Lx = " << Lx << "\n";
    cout << "t = " << t << "\n";
    cout << "U  = " << U << "\n";
    cout << "V  = " << V << "\n";
    if (params.ImpurityMode) {
      cout << "Impurity = " << params.Impurity << "\n";
    }
    cout << "total time = " << params.tau * params.steps << std::endl;

    cout << " ****** TDVP Parameter List ******" << "\n";
    cout << "tau = " << params.tau << "\n";
    cout << "steps = " << params.steps << "\n";
    cout << "Dmax = " << params.Dmax << "\n";
    cout << "CutOff = " << params.CutOff << "\n";
    cout << "LanczErr = " << params.LanczErr << std::endl;
  }
  clock_t startTime, endTime;
  startTime = clock();

  using namespace plain_hubbard;
  OperatorInitial();
  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  std::ifstream ifs("e0");
  double e0(0.0);
  if (ifs.good()) {
    ifs >> e0;
  } else {
    std::cout << "Open e0 file error. Please check the e0 file." << std::endl;
    exit(1);
  }
  ifs.close();
  std::cout << "e0 = " << e0 << std::endl;

  qlmps::TDVPEvolveParams<U1U1QN> sweep_params;
  std::string measure_file_base_name;

  Tensor op1, op2;
  switch (params.CorrelationMode) {
    case 0:op1 = sz;
      op2 = sz;
      measure_file_base_name = "szsz_dynamic";
      break;
    case 1:op1 = sp;
      op2 = sm;
      measure_file_base_name = "spsm_dynamic";
      break;
    case 2:op1 = sm;
      op2 = sp;
      measure_file_base_name = "smsp_dynamic";
      break;
    case 3:op1 = id;
      op2 = sz;
      measure_file_base_name = "szszc_dynamic"; //continue
      break;
    case 4:op1 = id;
      op2 = sm;
      measure_file_base_name = "spsmc_dynamic"; //continue
      break;
    case 5:op1 = id;
      op2 = sp;
      measure_file_base_name = "smspc_dynamic"; //continue
      break;
    default:std::cout << "Not support correlation mode. exit(1)" << std::endl;
      exit(1);
  }
  sweep_params = qlmps::TDVPEvolveParams(
      params.tau, params.steps,
      N / 2,
      op1, id, op2, id,
      e0,
      params.Dmin, params.Dmax, params.CutOff,
      qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter)
  );

  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_outF);
  qlmps::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);

  qlmps::MPO<Tensor> mpo(N);
  if (IsPathExist(kMpoPath)) {
    for (size_t i = 0; i < mpo.size(); i++) {
      std::string filename = kMpoPath + "/" +
          kMpoTenBaseName + std::to_string(i) + "." + kQLTenFileSuffix;
      mpo.LoadTen(i, filename);
    }

    cout << "MPO loaded." << endl;
  } else {
    cout << "No mpo directory. exiting" << std::endl;
    exit(0);
  }

  using FiniteMPST = qlmps::FiniteMPS<TenElemT, U1U1QN>;
  FiniteMPST mps(sites);

  if (IsPathExist(kMpsPath)) {
    if (N == GetNumofMps()) {
      cout << "The number of mps files is consistent with mps size." << endl;
      cout << "Directly use mps from files." << endl;
    } else {
      cout << "Mps file number is wrong." << endl;
      exit(1);
    }
  } else {
    cout << "No mps directory." << endl;
    exit(1);
  }

  if (world.size() == 1) {
    TwoSiteFiniteTDVP(mps, mpo, sweep_params, measure_file_base_name);
  } else {
    qlmps::MPITDVPSweepParams<U1U1QN> sweep_params;

    sweep_params = qlmps::MPITDVPSweepParams(
        params.tau, params.steps,
        N / 2,
        op1, id, op2, id,
        e0,
        params.Dmin, params.Dmax, params.CutOff,
        qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter)
    );

    TwoSiteFiniteTDVP(mps, mpo, sweep_params, measure_file_base_name, world);
  }

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
