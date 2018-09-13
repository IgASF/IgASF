/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <maps/transform_coefs.h>
#include <tools/timing.h>

namespace{
// implementation forward

template<int dom, int tar>
void transform(
               FMatrix &jacs,
               view<const real_t> A,
               view<const real_t> b,
               view<const real_t> c,
               std::vector<view<real_t>> &transformedA,
               std::vector<view<real_t>> &transformedB,
               view<real_t>              &transformedC
               );
} // anonymous namespace

// public interface

void transformCoefs(
                    const GeoPtr& geo,
                    const CartesianGrid &nodes,
                    view<const real_t> A,
                    view<const real_t> b,
                    view<const real_t> c,
                    std::vector<view<real_t>> tA,
                    std::vector<view<real_t>> tB,
                    view<real_t>              tC
                    )
{
    const int dom=geo->domainDim();
    const int tar=geo->targetDim();
    
    FMatrix jacs=geo->jacobian(nodes);
    
    auto t_start = std::chrono::high_resolution_clock::now();
    switch ( (dom-1)*4+(tar-1) )
    {
        // dom=1
        case 0: transform<1,1>( jacs, A, b, c, tA, tB,tC ); break;
        case 1: transform<1,2>( jacs, A, b, c, tA, tB,tC ); break;
        case 2: transform<1,3>( jacs, A, b, c, tA, tB,tC ); break;
        case 3: transform<1,4>( jacs, A, b, c, tA, tB,tC ); break;
        // dom=2
        case 5: transform<2,2>( jacs, A, b, c, tA, tB,tC ); break;
        case 6: transform<2,3>( jacs, A, b, c, tA, tB,tC ); break;
        case 7: transform<2,4>( jacs, A, b, c, tA, tB,tC ); break;
        // dom=3
        case 10: transform<3,3>( jacs, A, b, c, tA, tB,tC ); break;
        case 11: transform<3,4>( jacs, A, b, c, tA, tB,tC ); break;
        // dom 4
        case 15: transform<4,4>( jacs, A, b, c, tA, tB,tC ); break;
        default:
        assert(false);
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    time_geo_transform.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
}


void transformCoefs(
                    const GeoPtr& geo,
                    const CartesianGrid &nodes,
                    EqCoef eq,
                    std::vector<view<real_t>> transformedA,
                    std::vector<view<real_t>> transformedB,
                    view<real_t>              transformedC)
{
    int tar=geo->targetDim();
    int tar2=tar*tar;
    int num=nodes.numPoints();
    
    OwningView<real_t> Adata( eq.hasA() ? num*tar2: 0);
    OwningView<real_t> Bdata( eq.hasB() ? num*tar: 0);
    std::vector<real_t> Cdata( eq.hasC() ? num: 0, eq.C);
    
    view<real_t> Av(Adata);
    view<real_t> Bv(Bdata);
    view<real_t> Cv(Cdata);
    
    if (eq.hasA()) Av.matrix(tar2,num).colwise() = view<real_t>(eq.A).vector();
    if (eq.hasB()) Bv.matrix(tar,num).colwise()  = view<real_t>(eq.B).vector();
    
    transformCoefs(geo, nodes, Av,Bv,Cv,transformedA, transformedB, transformedC );
}

namespace {
// implementation

template <int dom, int tar>
void transformA( view<const real_t> jac, view<const real_t>  det, view<const real_t> A, std::vector<view<real_t>> &tA )
{
    const index_t blockJ=dom*tar;
    const index_t blockA=tar*tar;
    const index_t num=det.size();
    
    Eigen::Matrix<real_t,dom,dom> tmp;
    
    for (index_t p=0; p<num;++p)
        {
        auto Jp=jac.subView(blockJ*p,blockJ).matrix<dom,tar>();
        auto Ap=A.subView(p*blockA,blockA).matrix<tar,tar>();
        tmp=Jp*Ap*Jp.transpose()*det[p];
        //if (std::getenv("debug")) { std::cout<<"J^-t A J-1 det \n"<<tmp<<"\n"<<std::endl; }
        for (int i=0;i<dom;++i)
            for (int j=0;j<dom;++j)
                tA[i*dom+j][p]=tmp(i,j);
        }
}

template <int dom, int tar>
void transformB( view<const real_t> jac, view<const real_t>  det, view<const real_t> B, std::vector<view<real_t>> &tB )
{
    const index_t blockJ=dom*tar;
    const index_t num=det.size();
    
    Eigen::Matrix<real_t,dom,1> tmp;
    for (index_t p=0; p<num;++p)
        {
        auto Jp=jac.subView(p*blockJ,blockJ).matrix<dom,tar>();
        auto Bp=B.subView(p*tar,tar).matrix<tar,1>();
        tmp=Jp*Bp*det[p];
        for (int i=0;i<dom;++i)
            tB[i][p]=tmp(i);
        }
}

template <int dom, int tar>
void transformC(const view< real_t>  det,  view< const real_t> C, view<real_t> tC )
{
    const index_t num=det.size();
    tC.matrix(num,1).array()=det.matrix(num,1).array()*C.matrix(num,1).array();
}



template <int dom, int tar>
real_t invertInPlace (view<real_t> jac)
{
    real_t det;
    auto J =jac.matrix<tar,dom>();
    auto JI=jac.matrix<dom,tar>();
    
    if constexpr(dom==tar)
        {
        det=J.determinant();
        JI=J.inverse().eval();
        }
    else
        {
        Eigen::Matrix<real_t,dom,dom> tmp;
        tmp = J.transpose()*J;
        det=std::sqrt(tmp.determinant());
        JI=(J*tmp.inverse()).transpose();
        }
    return det;
}

template<int dom, int tar>
void transform(
               FMatrix &jacs,
               view<const real_t> A,
               view<const real_t> B,
               view<const real_t> C,
               std::vector<view<real_t>> &tA,
               std::vector<view<real_t>> &tB,
               view<real_t>              &tC
               )
{
    OwningView<real_t> det(jacs.cols());
    view<real_t> jacV(jacs);
    for (index_t p=0;p<jacs.cols();++p)
    {
        det[p]=invertInPlace<dom,tar>(jacV.subView(dom*tar*p,dom*tar));
    }
    
    if (A.size()) transformA<dom,tar>( jacV, det, A,  tA );
    if (B.size()) transformB<dom,tar>( jacV, det, B,  tB );
    if (C.size()) transformC<dom,tar>( det, C,  tC );
}
} // anonymous namespace
