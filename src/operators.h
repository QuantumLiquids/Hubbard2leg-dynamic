/*
 * File Name: operators.h
 * Description: Declare the operators in Hubbard DMRG project
 * Created by Hao-Xin on 2021/10/15.
 *
 */

#ifndef HUBBARD_SRC_OPERATORS_H_
#define HUBBARD_SRC_OPERATORS_H_

#include "qldouble.h"

namespace plain_hubbard {
//Fermionic operators
extern Tensor sz, sp,sm, id;
extern Tensor f, bupc, bupa, bdnc, bdna;
extern Tensor bupcF, bupaF, Fbdnc, Fbdna, bdncF;
extern Tensor cupccdnc, cdnacupa, Uterm, nf, nfsquare, nup,ndn;
void OperatorInitial();
}

namespace z4_hubbard {
extern std::vector<Tensor> sz, sp,sm, id,
    f, bupc, bupa, bdnc, bdna,
    bupcF, bupaF, Fbdnc, Fbdna,
    bdncF,
    Uterm,
    n,
    cupccdnc, cdnacupa;
void OperatorInitial();
}

namespace z6_hubbard {
extern std::vector<Tensor> sz, sp,sm, id,
    f, bupc, bupa, bdnc, bdna,
    bupcF, bupaF, Fbdnc, Fbdna,
    bdncF,
    Uterm,
    n,
    cupccdnc, cdnacupa;
void OperatorInitial();
}



#endif //HUBBARD_SRC_OPERATORS_H_
