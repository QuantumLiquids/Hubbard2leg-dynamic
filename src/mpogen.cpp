/*
 * File Name: mpogen.cpp
 * Description: Generate and dump mpo for 1D Hubbard model
 * Created by Hao-Xin on 2021/11/07.
 *
 */

#include "qldouble.h"
#include "operators.h"
#include "params_case.h"
#include "qlmps/qlmps.h"

#include <time.h>

using namespace std;
using namespace qlten;
using namespace qlmps;

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);
  const size_t Lx = params.Lx;
  const size_t N = 2 * Lx;
  cout << "The total number of sites: " << N << endl;
  double t = params.t, U = params.U, V = params.V;
  double t2 = params.t2;
  double t_perp = params.t_perp;
  using std::cout;
  cout << " ****** 1D Hubbard Model Parameter List ****** " << "\n";
  cout << "Lx = " << Lx << "\n";
  cout << "t = " << t << "\n";
  cout << "tp = " << t_perp << "\n";
  cout << "t2 = " << t2 << "\n";
  cout << "U  = " << U << "\n";
  cout << "V  = " << V << "\n";

  clock_t startTime, endTime;
  startTime = clock();

  if (!IsPathExist(kMpoPath)) {
    CreatPath(kMpoPath);
  }

  using namespace plain_hubbard;
  OperatorInitial();
  const SiteVec<TenElemT, U1U1QN> sites = SiteVec<TenElemT, U1U1QN>(N, pb_outF);
  qlmps::MPOGenerator<TenElemT, U1U1QN> mpo_gen(sites, qn0);

  for (size_t i = 0; i < N; ++i) {
    mpo_gen.AddTerm(U, Uterm, i);
    cout << "add site" << i << "Hubbard U term" << endl;
  }

  //NN hopping terms and density-density interaction
  for (size_t i = 0; i < N - 3; i = i + 2) {
    //horizontal hopping and attraction
    size_t site1 = i, site2 = i + 2;
    mpo_gen.AddTerm(-t, bupcF, site1, bupa, site2, f);
    mpo_gen.AddTerm(-t, bdnc, site1, Fbdna, site2, f);
    mpo_gen.AddTerm(t, bupaF, site1, bupc, site2, f);
    mpo_gen.AddTerm(t, bdna, site1, Fbdnc, site2, f);
    mpo_gen.AddTerm(V, nf, site1, nf, site2);
    cout << "add site (" << site1 << "," << site2 << ")  hopping and density interaction terms" << endl;

    site1 = i + 1, site2 = i + 3;
    mpo_gen.AddTerm(-t, bupcF, site1, bupa, site2, f);
    mpo_gen.AddTerm(-t, bdnc, site1, Fbdna, site2, f);
    mpo_gen.AddTerm(t, bupaF, site1, bupc, site2, f);
    mpo_gen.AddTerm(t, bdna, site1, Fbdnc, site2, f);
    mpo_gen.AddTerm(V, nf, site1, nf, site2);
    cout << "add site (" << site1 << "," << site2 << ")  hopping and density interaction terms" << endl;

    //diagonal hopping
    site1 = i, site2 = i + 3;
    mpo_gen.AddTerm(-t2, bupcF, site1, bupa, site2, f);
    mpo_gen.AddTerm(-t2, bdnc, site1, Fbdna, site2, f);
    mpo_gen.AddTerm(t2, bupaF, site1, bupc, site2, f);
    mpo_gen.AddTerm(t2, bdna, site1, Fbdnc, site2, f);

    site1 = i + 1, site2 = i + 2;
    mpo_gen.AddTerm(-t2, bupcF, site1, bupa, site2, f);
    mpo_gen.AddTerm(-t2, bdnc, site1, Fbdna, site2, f);
    mpo_gen.AddTerm(t2, bupaF, site1, bupc, site2, f);
    mpo_gen.AddTerm(t2, bdna, site1, Fbdnc, site2, f);

//    mpo_gen.AddTerm(J, sp, site1, sm, site2);
//    mpo_gen.AddTerm(J, sm, site1, sp, site2);
//    mpo_gen.AddTerm(2 * J, sz, site1, sz, site2);
//
//    mpo_gen.AddTerm(-J, cupccdnc, site1, cdnacupa, site2);
//    mpo_gen.AddTerm(-J, cdnacupa, site1, cupccdnc, site2);
//    mpo_gen.AddTerm(J / 2, nf + (-id), site1, nf + (-id), site2);
//    cout << "add site (" << site1 << "," << site2 << ")  hopping square terms" << endl;
  }

  //vertical hopping and attraction
  for (size_t i = 0; i < N; i = i + 2) {
    size_t site1 = i, site2 = i + 1;
    mpo_gen.AddTerm(-t_perp, bupcF, site1, bupa, site2, f);
    mpo_gen.AddTerm(-t_perp, bdnc, site1, Fbdna, site2, f);
    mpo_gen.AddTerm(t_perp, bupaF, site1, bupc, site2, f);
    mpo_gen.AddTerm(t_perp, bdna, site1, Fbdnc, site2, f);
    mpo_gen.AddTerm(V, nf, site1, nf, site2);
    cout << "add site (" << site1 << "," << site2 << ")  hopping and density interaction terms" << endl;
  }
  const size_t tdvp_reference_point = N / 2;
  // assume N%4 == 0
  // params.Impurity > 0, means trapping hole on the site.
  switch (params.ImpurityMode) {
    case 0: {
      break;
    }
    case 1 : {
      mpo_gen.AddTerm(params.Impurity, nf, N / 2 + 1);
      break;
    }
    case 2: {
      mpo_gen.AddTerm(params.Impurity, nf, N / 2 - 2);
      break;
    }
    case 3: {
      mpo_gen.AddTerm(params.Impurity, nf, N / 2 - 1);
      break;
    }
  }
//
//  double Perturbation = params.Perturbation;
//  if (fabs(Perturbation) > kDoubleEpsilon) {
//    size_t ChargePeriod = params.wavelength;//suppose be 4
//    if (L % params.wavelength == 0) {
//      for (size_t i = 0; i < L; i++) {
//        size_t x = i;
//        double amplitude = -Perturbation * cos(M_PI / ChargePeriod + x * (2 * M_PI / ChargePeriod));
//        mpo_gen.AddTerm(amplitude, nf, i);
//      }
//    } else if (L % params.wavelength == 1) { //Lx = 5, 9, 11
//      for (size_t i = 0; i < L; i++) {
//        size_t x = i;
//        double amplitude = -Perturbation * cos(x * (2 * M_PI / ChargePeriod));
//        mpo_gen.AddTerm(amplitude, nf, i);
//      }
//    } else {
//      std::cout << "L setting is not support for add charge density perturbation!" << std::endl;
//      exit(1);
//    }
//
//    cout << "Add charge density pinning field, mu = " << Perturbation << ", with Period = " << ChargePeriod << endl;
//  }

  auto mpo = mpo_gen.Gen();
  cout << "MPO generated." << endl;

  for (size_t i = 0; i < mpo.size(); i++) {
    std::string filename = kMpoPath + "/" +
        kMpoTenBaseName + std::to_string(i)
        + "." + kQLTenFileSuffix;
    mpo.DumpTen(i, filename);
  }

  endTime = clock();
  cout << "CPU Time : " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}




