/*
 * File Name: operators.cpp
 * Description: Define the operators in Hubbard DMRG project
 * Created by Hao-Xin on 2021/10/15.
 *
 */

#include "qldouble.h"

namespace plain_hubbard {
Tensor sz = Tensor({pb_inF, pb_outF});
Tensor sp = Tensor({pb_inF, pb_outF});
Tensor sm = Tensor({pb_inF, pb_outF});
Tensor id = Tensor({pb_inF, pb_outF});

Tensor f = Tensor({pb_inF,pb_outF}); //fermion's insertion operator

auto bupc = Tensor({pb_inF,pb_outF}); //hardcore boson, b_up^creation, used for JW transformation
auto bupa = Tensor({pb_inF,pb_outF}); //hardcore boson, b_up^annihilation
auto bdnc = Tensor({pb_inF,pb_outF}); //hardcore boson, b_down^creation
auto bdna = Tensor({pb_inF,pb_outF}); //hardcore boson, b_down^annihilation


auto bupcF = Tensor({pb_inF,pb_outF}); // matrix product of bupc * f
auto bupaF = Tensor({pb_inF,pb_outF});
auto Fbdnc = Tensor({pb_inF,pb_outF});
auto Fbdna = Tensor({pb_inF,pb_outF});
auto bdncF = Tensor({pb_inF,pb_outF});


auto cupccdnc = Tensor({pb_inF,pb_outF}); // c_up^creation * c_down^creation=b_up^creation*b_down^creation*F

auto cdnacupa = Tensor({pb_inF,pb_outF}); // onsite pair, usually c_up*c_dn


auto Uterm = Tensor({pb_inF,pb_outF}); // Hubbard Uterm, nup*ndown

auto nf = Tensor({pb_inF,pb_outF}); // nup+ndown, fermion number

auto nfsquare = Tensor({pb_inF,pb_outF}); // nf^2
auto nup = Tensor({pb_inF,pb_outF}); // fermion number of spin up
auto ndn = Tensor({pb_inF,pb_outF}); // ndown


void OperatorInitial(){
  static bool initialized = false;
  if(!initialized){
    /* Since the first indexes of operators will be contracted to the wave function |Psi>,
     * so we should regard the first indexes as the row index
     *
     */

    sz({1, 1}) = 0.5;
    sz({2, 2}) = -0.5;
    sp({2, 1}) = 1.0;
    sm({1, 2}) = 1.0;
    id({0,0}) = 1;
    id({1,1}) = 1;
    id({2,2}) = 1;
    id({3,3}) = 1;

    f({0,0}) = 1;
    f({1,1}) = -1;
    f({2,2}) = -1;
    f({3,3}) = 1;


    bupc({2,0}) = 1;
    bupc({3,1}) = 1;
    bdnc({1,0}) = 1;
    bdnc({3,2}) = 1;
    bupa({0,2}) = 1;
    bupa({1,3}) = 1;
    bdna({0,1}) = 1;
    bdna({2,3}) = 1;


    bupcF({2,0}) = -1;
    bupcF({3,1}) = 1;
    Fbdnc({1,0}) = 1;
    Fbdnc({3,2}) = -1;
    bupaF({0,2}) = 1;
    bupaF({1,3}) = -1;
    Fbdna({0,1}) = -1;
    Fbdna({2,3}) = 1;

    bdncF = -Fbdnc;

    cupccdnc({3,0}) = 1;
    cdnacupa({0,3}) = 1;


    Uterm({0,0}) = 1;


    nf({0,0}) = 2;
    nf({1,1}) = 1;
    nf({2,2}) = 1;


    nfsquare({0,0}) = 4;
    nfsquare({1,1}) = 1;
    nfsquare({2,2}) = 1;

    nup({0,0}) = 1;
    nup({1,1}) = 1;
    ndn({0,0}) = 1;
    ndn({2,2}) = 1;

    initialized=true;
  }
}

}//namespace plain_hubbard


namespace z4_hubbard {
using OpSet = std::vector<Tensor>;
OpSet sz = OpSet(Ly);
OpSet sp = sz, sm = sz, id = sz;
OpSet f = sz;
OpSet bupc=sz, bupa=sz, bdnc=sz, bdna=sz;
OpSet bupcF=sz, bupaF=sz, Fbdnc=sz, Fbdna=sz;
OpSet bdncF=sz;
OpSet Uterm = sz;
OpSet n = sz;
OpSet cupccdnc = sz, cdnacupa = sz;


void OperatorInitial() {
  /*
   * Here the first indexes of operator remains to corresponding to the column indexes, should be changed if using in future.
   */
  static bool initialized = false;
  if(initialized) return;//this function only run one time

  if(!hilbert_space_species_initialed){
    InitialHilbertSpaceSpecies();
  }
  for(size_t i = 0; i < Ly; i++){
    sz[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    sp[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    sm[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    id[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});

    sz[i]({1, 1}) = 0.5;
    sz[i]({2, 2}) = -0.5;
    sp[i]({1, 2}) = 1.0;
    sm[i]({2, 1}) = 1.0;
    id[i]({0,0}) = 1;
    id[i]({1,1}) = 1;
    id[i]({2,2}) = 1;
    id[i]({3,3}) = 1;

    f[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    f[i]({0,0}) = 1;
    f[i]({1,1}) = -1;
    f[i]({2,2}) = -1;
    f[i]({3,3}) = 1;

    bupc[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bupa[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bdnc[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bdna[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});

    bupc[i]({0,2}) = 1;
    bupc[i]({1,3}) = 1;
    bdnc[i]({0,1}) = 1;
    bdnc[i]({2,3}) = 1;
    bupa[i]({2,0}) = 1;
    bupa[i]({3,1}) = 1;
    bdna[i]({1,0}) = 1;
    bdna[i]({3,2}) = 1;

    bupcF[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bupaF[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    Fbdnc[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    Fbdna[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});

    bupcF[i]({0,2}) = -1;
    bupcF[i]({1,3}) = 1;
    Fbdnc[i]({0,1}) = 1;
    Fbdnc[i]({2,3}) = -1;
    bupaF[i]({2,0}) = 1;
    bupaF[i]({3,1}) = -1;
    Fbdna[i]({1,0}) = -1;
    Fbdna[i]({3,2}) = 1;

    bdncF[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bdncF[i]({0,1}) = -1;
    bdncF[i]({2,3}) = 1;

    Uterm[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});// nup*ndown
    Uterm[i]({0,0}) = 1;

    n[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});// nup+ndown
    n[i]({0,0}) = 2;
    n[i]({1,1}) = 1;
    n[i]({2,2}) = 1;


    cupccdnc[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});// c_up^creation * c_down^creation=b_up^creation*b_down^creation*F
    cupccdnc[i]({0,3}) = 1;

    cdnacupa[i] = Tensor({pb_inF_species[i], pb_outF_species[i]}); // onsite pair, usually c_up*c_dn
    cdnacupa[i]({3,0}) = 1;
  }
  initialized=true;
}
}//z4_hubbard


namespace  z6_hubbard {
//completely copy from z4_hubbard
using OpSet = std::vector<Tensor>;
OpSet sz = OpSet(Ly);
OpSet sp = sz, sm = sz, id = sz;
OpSet f = sz;
OpSet bupc=sz, bupa=sz, bdnc=sz, bdna=sz;
OpSet bupcF=sz, bupaF=sz, Fbdnc=sz, Fbdna=sz;
OpSet bdncF=sz;
OpSet Uterm = sz;
OpSet n = sz;
OpSet cupccdnc = sz, cdnacupa = sz;

void OperatorInitial() {
  static bool initialized = false;
  if(initialized) return;//this function only run one time

  if(!hilbert_space_species_initialed){
    InitialHilbertSpaceSpecies();
  }
  for(size_t i = 0; i < Ly; i++){
    sz[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    sp[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    sm[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    id[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});

    sz[i]({1, 1}) = 0.5;
    sz[i]({2, 2}) = -0.5;
    sp[i]({1, 2}) = 1.0;
    sm[i]({2, 1}) = 1.0;
    id[i]({0,0}) = 1;
    id[i]({1,1}) = 1;
    id[i]({2,2}) = 1;
    id[i]({3,3}) = 1;

    f[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    f[i]({0,0}) = 1;
    f[i]({1,1}) = -1;
    f[i]({2,2}) = -1;
    f[i]({3,3}) = 1;

    bupc[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bupa[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bdnc[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bdna[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});

    bupc[i]({0,2}) = 1;
    bupc[i]({1,3}) = 1;
    bdnc[i]({0,1}) = 1;
    bdnc[i]({2,3}) = 1;
    bupa[i]({2,0}) = 1;
    bupa[i]({3,1}) = 1;
    bdna[i]({1,0}) = 1;
    bdna[i]({3,2}) = 1;

    bupcF[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bupaF[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    Fbdnc[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    Fbdna[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});

    bupcF[i]({0,2}) = -1;
    bupcF[i]({1,3}) = 1;
    Fbdnc[i]({0,1}) = 1;
    Fbdnc[i]({2,3}) = -1;
    bupaF[i]({2,0}) = 1;
    bupaF[i]({3,1}) = -1;
    Fbdna[i]({1,0}) = -1;
    Fbdna[i]({3,2}) = 1;

    bdncF[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});
    bdncF[i]({0,1}) = -1;
    bdncF[i]({2,3}) = 1;

    Uterm[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});// nup*ndown
    Uterm[i]({0,0}) = 1;

    n[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});// nup+ndown
    n[i]({0,0}) = 2;
    n[i]({1,1}) = 1;
    n[i]({2,2}) = 1;


    cupccdnc[i] = Tensor({pb_inF_species[i], pb_outF_species[i]});// c_up^creation * c_down^creation=b_up^creation*b_down^creation*F
    cupccdnc[i]({0,3}) = 1;

    cdnacupa[i] = Tensor({pb_inF_species[i], pb_outF_species[i]}); // onsite pair, usually c_up*c_dn
    cdnacupa[i]({3,0}) = 1;
  }
  initialized=true;
}

}//z6_hubbard






