/*
 * File Name: dmrg.cpp
 * Description: Calculate the ground state for 1D Hubbard model
 * Created by Hao-Xin on 2021/11/07.
 *
 * Usage:
 *    mpirun -np <the number of process> ./dmrg params.json
 * Optional parameter:
 *    --D=100,200,5000...., will cover the max bond dimension in params.json
 */


#include <iostream>
#include "qlmps/qlmps.h"
#include "qlten/qlten.h"
#include <time.h>
#include "params_case.h"
#include "qldouble.h"
#include "myutil.h"

using namespace qlmps;
using namespace qlten;
using namespace std;

bool ParserBondDimension(int argc, char *argv[],
                         vector<size_t> &D_set);

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env;
  mpi::communicator world;
  std::vector<size_t> output_D_set;
  bool has_bond_dimension_parameter = ParserBondDimension(
      argc, argv,
      output_D_set);

  CaseParams params(argv[1]);
  size_t Lx = params.Lx;
  size_t N = 2 * Lx;

  double t = params.t, U = params.U, V = params.V;
  using std::cout;
  cout << " ****** 1D Hubbard Model Parameter List ****** " << "\n";
  cout << "Lx = " << Lx << "\n";
  cout << "t = " << t << "\n";
  cout << "U  = " << U << "\n";
  cout << "V  = " << V << "\n";
  if (params.ImpurityMode) {
    cout << "Impurity  = " << params.Impurity << "\n";
  }
  cout << " ****** DMRG Parameter List ******" << "\n";
  cout << "Sweeps = " << params.Sweeps << "\n";
  cout << "Dmax = " << params.Dmax << "\n";
  cout << "CutOff = " << params.CutOff << "\n";
  cout << "LanczErr = " << params.LanczErr << std::endl;

  clock_t startTime, endTime;
  startTime = clock();

  qlten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  bool noise_valid(false);
  for (size_t i = 0; i < params.noise.size(); i++) {
    if (params.noise[i] != 0) {
      noise_valid = true;
      break;
    }
  }

  double e0(0.0); //energy


  using namespace plain_hubbard;
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

  std::vector<long unsigned int> stat_labs(N);
  if (params.Numhole == 0) {
    for (size_t i = 0; i < N; i++) {
      stat_labs[i] = i % 2 + 1;
    }
  } else if (N % params.Numhole == 0) {
    const size_t site_number_per_hole = N / params.Numhole;
    int sz_label = 0;
    for (size_t i = 0; i < N; ++i) {
      if (i % site_number_per_hole == site_number_per_hole - 1) {
        stat_labs[i] = 3;
      } else {
        stat_labs[i] = sz_label % 2 + 1;
        sz_label++;
      }
    }
  } else {
    for (size_t i = 0; i < params.Numhole; i++) {
      stat_labs[i] = 3;
    }
    for (size_t i = params.Numhole; i < N; i++) {
      stat_labs[i] = i % 2 + 1;
    }
  }

  qlmps::FiniteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter)
  );

  if (world.rank() == 0) {
    if (IsPathExist(kMpsPath)) {
      if (N == GetNumofMps()) {
        cout << "The number of mps files is consistent with mps size." << endl;
        cout << "Directly use mps from files." << endl;
      } else {
        qlmps::DirectStateInitMps(mps, stat_labs);
        cout << "Initial mps as direct product state." << endl;
        mps.Dump(sweep_params.mps_path, true);
      }
    } else {
      qlmps::DirectStateInitMps(mps, stat_labs);
      cout << "Initial mps as direct product state." << endl;
      mps.Dump(sweep_params.mps_path, true);
    }
  }

  if (!has_bond_dimension_parameter) {
    if (world.size() == 1) {
      e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
    } else {
      qlmps::FiniteVMPSSweepParams sweep_params(
          params.Sweeps,
          params.Dmin, params.Dmax, params.CutOff,
          qlmps::LanczosParams(params.LanczErr, params.MaxLanczIter)
      );
      e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
    }
  } else {
    size_t DMRG_time = output_D_set.size();
    std::vector<size_t> MaxLanczIterSet(DMRG_time);
    MaxLanczIterSet.back() = params.MaxLanczIter;
    if (DMRG_time > 1) {
      size_t MaxLanczIterSetSpace;
      MaxLanczIterSet[0] = 3;
      MaxLanczIterSetSpace = (params.MaxLanczIter - 3) / (DMRG_time - 1);
      std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << ", ";
      for (size_t i = 1; i < DMRG_time - 1; i++) {
        MaxLanczIterSet[i] = MaxLanczIterSet[i - 1] + MaxLanczIterSetSpace;
        std::cout << MaxLanczIterSet[i] << ", ";
      }
      std::cout << MaxLanczIterSet.back() << "]" << std::endl;
    } else {
      std::cout << "Setting MaxLanczIter as : [" << MaxLanczIterSet[0] << "]" << std::endl;
    }

    if (world.size() == 1) {
      for (size_t i = 0; i < DMRG_time; i++) {
        size_t D = output_D_set[i];
        std::cout << "D_max = " << D << std::endl;
        qlmps::FiniteVMPSSweepParams sweep_params(
            params.Sweeps,
            D, D, params.CutOff,
            qlmps::LanczosParams(params.LanczErr, MaxLanczIterSet[i])
        );

        e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
      }
    } else {
      for (size_t i = 0; i < DMRG_time; i++) {
        size_t D = output_D_set[i];
        std::cout << "D_max = " << D << std::endl;
        qlmps::FiniteVMPSSweepParams sweep_params(
            params.Sweeps,
            D, D, params.CutOff,
            qlmps::LanczosParams(params.LanczErr, MaxLanczIterSet[i])
        );
        e0 = qlmps::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
      }
    }
  }

  std::cout << "E0/site: " << e0 / N << std::endl;

  if (world.rank() == 0) {
    std::ofstream file("e0");
    file << e0 << std::endl;
    file.close();
  }

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}

