/*
 * File Name: gqdouble.h
 * Description: Define the types and constant in Hubbard DMRG project
 * Created by Hao-Xin on 2021/10/15.
 *
 */

#ifndef HUBBARD_SRC_GQDOUBLE_H_
#define HUBBARD_SRC_GQDOUBLE_H_

#include "qlten/qlten.h"


using TenElemT = qlten::QLTEN_Complex;


using U1U1QN = qlten::special_qn::U1U1QN;
using U1U1Z4QN = qlten::special_qn::U1U1ZnQN<4>;
using U1U1Z6QN = qlten::special_qn::U1U1ZnQN<6>;
using U1U1Z8QN = qlten::special_qn::U1U1ZnQN<8>;

using qlten::QLTensor;


namespace plain_hubbard {
using QNSctT = qlten::QNSector<U1U1QN>;
using IndexT = qlten::Index<U1U1QN>;
using Tensor = QLTensor<TenElemT, U1U1QN>;
const U1U1QN qn0 = U1U1QN(0,0); //N(particle number), Sz
const IndexT pb_outF = IndexT({   //QNSctT( U1U1QN(N, Sz), degeneracy )
                                  QNSctT(U1U1QN(2,0), 1),
                                  QNSctT(U1U1QN(1,1), 1),
                                  QNSctT(U1U1QN(1, -1), 1),
                                  QNSctT(U1U1QN(0,0), 1)
                              },
                              qlten::TenIndexDirType::OUT
);
const IndexT pb_inF = qlten::InverseIndex(pb_outF);

}//plain_hubbard

namespace z4_hubbard {
const size_t Ly = 4;
using QNSctT = qlten::QNSector<U1U1Z4QN>;
using IndexT = qlten::Index<U1U1Z4QN>;
using IndexVecT = qlten::IndexVec<U1U1Z4QN >;
using Tensor = QLTensor<TenElemT, U1U1Z4QN>;
const U1U1Z4QN qn0 = U1U1Z4QN(0, 0, 0); //N, Sz, k
inline IndexVecT pb_outF_species = IndexVecT(Ly);
inline IndexVecT pb_inF_species = IndexVecT(Ly);

inline bool hilbert_space_species_initialed = false;
inline void InitialHilbertSpaceSpecies(){
  if(hilbert_space_species_initialed) return;
  for (size_t i = 0; i < Ly; i++){
    //i is the coordinate in y direction(momentum space actually)
    pb_outF_species[i] = IndexT({
                                    QNSctT(U1U1Z4QN(2, 0, (2*i)%Ly), 1),
                                    QNSctT(U1U1Z4QN(1, 1, i), 1),
                                    QNSctT(U1U1Z4QN(1, -1, i), 1),
                                    QNSctT(U1U1Z4QN(0, 0, 0), 1)
                                },
                                qlten::TenIndexDirType::OUT
    );
    pb_inF_species[i] = InverseIndex( pb_outF_species[i] );
  }
  hilbert_space_species_initialed = true;
}

}//z4_hubbard

namespace z6_hubbard {
const size_t Ly = 6;
using QNSctT = qlten::QNSector<U1U1Z6QN>;
using IndexT = qlten::Index<U1U1Z6QN>;
using IndexVecT = qlten::IndexVec<U1U1Z6QN >;
using Tensor = QLTensor<TenElemT, U1U1Z6QN>;
const U1U1Z6QN qn0 = U1U1Z6QN(0, 0, 0); //N, Sz, k
inline IndexVecT pb_outF_species = IndexVecT(Ly);
inline IndexVecT pb_inF_species = IndexVecT(Ly);

inline bool hilbert_space_species_initialed = false;
inline void InitialHilbertSpaceSpecies(){
  if(hilbert_space_species_initialed) return;
  for (size_t i = 0; i < Ly; i++){
    //i is the coordinate in y direction(momentum space actually)
    pb_outF_species[i] = IndexT({
                                    QNSctT(U1U1Z6QN(2, 0, (2*i)%Ly), 1),
                                    QNSctT(U1U1Z6QN(1, 1, i), 1),
                                    QNSctT(U1U1Z6QN(1, -1, i), 1),
                                    QNSctT(U1U1Z6QN(0, 0, 0), 1)
                                },
                                qlten::TenIndexDirType::OUT
    );
    pb_inF_species[i] = InverseIndex( pb_outF_species[i] );
  }
  hilbert_space_species_initialed = true;
}

}//z6_hubbard





#endif //HUBBARD_SRC_GQDOUBLE_H_
