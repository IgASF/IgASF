/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <assembling/macroelement.h>
#include <bases/tensor_basis.h>
#include <bases/bspline.h>

#include <quadrature/per_element.h>
#include <quadrature/tensor_quadrature.h>

#include <tools/thread_pool.h>
#include <tools/multi_index.h>

#include <assembling/model.h>

#include <tools/timing.h>
#include <iostream>

namespace macroelement {
    
namespace {
    
    MMatrix initMatrix(const TensorBasis *test,const TensorBasis *trial, const TensorQuadrature  *quad)
    {
    auto dim=test->domainDim();
    std::vector<std::vector<PartialDerivative>> req(dim, std::vector<PartialDerivative>(1,PartialDerivative(0)) );
    // compute values so that we can compute the sparsity
    std::vector<BasisValues> tst=test->evaluateComponents (req,quad->grid());
    std::vector<BasisValues> trl=trial->evaluateComponents(req,quad->grid());
    // Precompute sparsity
    Sparsity krn, spr, tmp; // kroneker sparsity
    krn=bilinearSparsity(tst[0],trl[0]);
    for (int c=1; c<dim;++c)
        {
        std::swap(krn,tmp);
        spr = bilinearSparsity(tst[c],trl[c]);
        krn = kroneckerSparsity(spr,tmp);
        }
    return MMatrix(std::move(krn));
    }
    
    struct MacroInfo
    {
    Bspline tst;
    index_t tst_shift;
    Bspline trl;
    index_t  trl_shift;
    ElementQuadrature  quad;
    };
    typedef std::vector<MacroInfo>  Macro1D;
    typedef std::vector<Macro1D>    MacroPartition;
    
    Macro1D makeMacros(
                       const Bspline &test,
                       const Bspline &trial,
                       const ElementQuadrature &quad,
                       index_t macroS)
    {
    Macro1D result;
    
    // Compute element splitting
    auto &nodes   = quad.nodes();
    auto begTest  = ++(test.breaks().begin());
    auto endTest  = test.breaks().end();
    auto begTrial = ++( trial.breaks().begin());
    auto endTrial = trial.breaks().end();
    
    index_t pto=0;
    index_t pt=0;
    index_t ec=0;
    
    auto closeMacro=[&]()
        {
        // make quadrature
        view<const real_t> ns=view<const real_t>(quad.nodes()  ).subView(pto,pt-pto);
        view<const real_t> ws=view<const real_t>(quad.weights()).subView(pto,pt-pto);
        pto=pt;
        
        // make bases
        auto makeBasis =[&ns](const Bspline &bas, index_t beg)->Bspline
            {
            int d=bas.degree();
            index_t end = bas.evaluate({0}, ns.subView(ns.size()-1,1)).rowPos(0).back()+d+2; // index of last active + d+2
            view<const real_t> kns=view<const real_t>(bas.knots()).subView(beg,end-beg);
            return  Bspline(d,kns) ;
            };
        index_t tst_shift= test.evaluate ({0}, ns.subView(0,1)).rowPos(0).front(); // index of first active
        index_t trl_shift= trial.evaluate({0}, ns.subView(0,1)).rowPos(0).front(); // index of first active
        
        result.push_back(MacroInfo{
            makeBasis(test, tst_shift), tst_shift,
            makeBasis(trial,trl_shift), trl_shift,
            ElementQuadrature(ns,ws)
        });
        };
    
    while (begTest !=endTest && begTrial !=endTrial)
        {
        while( pt< quad.size() && nodes(0,pt) < *begTest &&  nodes(0,pt) < *begTrial)
            ++pt;
        ++ec;
        if (pt==quad.size() )  { closeMacro(); break;}
        if (ec==macroS)        { closeMacro(); ec=0; }
        while ( begTest  !=endTest  && nodes(0,pt) >= *begTest )  ++begTest;
        while ( begTrial !=endTrial && nodes(0,pt) >= *begTrial ) ++begTrial;
        }
    return result;
    }
    
    
    MacroPartition makeMacros ( const TensorBasis *test, const TensorBasis *trial, const TensorQuadrature  *quad, IMatrix macroSize)
    {
    const int dim = test->domainDim();
    MacroPartition result(dim);
    
    if (macroSize.rows()!=dim)
    {
        if (macroSize.rows()==0)
        {
            macroSize.resize(dim,1);
            macroSize.setConstant(-1);
        }
        else if (macroSize.rows()==1)
        {
            double tmp = macroSize(0,0);
            macroSize.resize(dim,1);
            macroSize.setConstant(tmp);
        }
        else
            throw std::runtime_error("The macro size does not agree with the dimension of the problem.");
    }
    
    for (int i=0;i<dim;++i )
    {
        if (macroSize(i)==-1)
            macroSize(i)= std::max((*test)[i].degree(), (*trial)[i].degree())+1;
        if (macroSize(i)<1)
            throw std::runtime_error("The macro size must be positive.");
    }

    FMatrix points(0,0);
    for (int i=0;i<dim;++i )
        {
        result[i] = makeMacros(
                               dynamic_cast<const Bspline&>((*test)[i]),
                               dynamic_cast<const Bspline&>((*trial)[i]),
                               *dynamic_cast<const ElementQuadrature*>((*quad)[i].get()), macroSize(i) );
        }
    return result;
    }
    
    
#define macroInfS 9
#define macroNum macroInf.col(0)
#define macroCom macroInf.col(1)
#define macroPos macroInf.col(2)
#define macroTST macroInf.col(3)
#define macroTRL macroInf.col(4)
#define macroTSTSUB macroInf.rightCols<4>().leftCols<2>()
#define macroTRLSUB macroInf.rightCols<2>()
    
    
    template <int DIM>
    void addPart(MMatrix &res, const MMatrix &part, const IMatrix &macroInf)
    {
    auto tstS=MultiIndex<DIM>(macroTST).subdomain(macroTSTSUB);
    auto trlS=MultiIndex<DIM>(macroTRL).subdomain(macroTRLSUB);
       
    for (index_t r=0;r<part.rows();++r)
        {
        index_t rr=tstS.local2global(r);
        
        const real_t  *srcD = part.data()+part.start(r);
        const auto     srcI = part.rowPos(r);
        
        real_t  *dstD = res.data()+res.start(rr);
        auto     dstI = res.rowPos(rr);
        
        for (index_t nz=0; nz<part.nnzs(r); ++nz )
            {
            index_t cc=trlS.local2global(srcI[nz]);
            while (dstI.front() < cc) { ++dstD;  dstI.popFront(); }
            *dstD+= srcD[nz];
            }
        }
    }
    
    
    void assembleMacroelement(const Model &eq, const MacroPartition &macros, MMatrix &result,  IMatrix macroInf)//, std::mutex &write)
    {
    const int domDim=macroInf.rows();
    
    std::vector<QuadPtr>  dirQ(domDim);
    std::vector<BasisPtr> dirBtst(domDim);
    std::vector<BasisPtr> dirBtrl(domDim);
    
    for (int i=0;i<domDim;++i)
        {
        dirQ[i]=QuadPtr(new ElementQuadrature(  macros[i][macroPos(i)].quad ));
        dirBtst[i]=BasisPtr(new Bspline( macros[i][macroPos(i)].tst ));
        dirBtrl[i]=BasisPtr(new Bspline( macros[i][macroPos(i)].trl ));
        macroTSTSUB(i,0) = macros[i][macroPos(i)].tst_shift;
        macroTSTSUB(i,1) = macroTSTSUB(i,0)+macros[i][macroPos(i)].tst.size();
        macroTRLSUB(i,0) = macros[i][macroPos(i)].trl_shift;
        macroTRLSUB(i,1) = macroTRLSUB(i,0)+macros[i][macroPos(i)].trl.size();
        }
    
    
    static thread_local MMatrix part;
    TensorBasis mtst(dirBtst);
    TensorBasis mtrl(dirBtrl);
    TensorQuadrature mquad(dirQ);
    
    eq.assemble(mtst, mtrl, mquad, part);

    DoNotOptimize(part);
    auto t_start = std::chrono::high_resolution_clock::now();
    switch (domDim)
        {
            case 4: addPart<4>(result,part,macroInf); break;
            case 3: addPart<3>(result,part,macroInf); break;
            case 2: addPart<2>(result,part,macroInf); break;
            case 1: addPart<1>(result,part,macroInf); break;
        }
    auto t_end = std::chrono::high_resolution_clock::now();
    time_add_macro.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    }

} // anonymous namespace
    
    MMatrix assembleParallel(const Model& eq, const BasisPtr &test, const BasisPtr &trial, int maxThread, IMatrix macroSize)
    {
    auto t_start = std::chrono::high_resolution_clock::now();

    const TensorBasis* TEST=dynamic_cast<const TensorBasis *>(test.get());
    const TensorBasis* TRIAL=dynamic_cast<const TensorBasis *>(trial.get());
    if (TEST==nullptr || TRIAL == nullptr  )
        throw std::logic_error("Not Implemented!");
    
    
    const int domDim = test->domainDim();
    QuadPtr quad = getRecommendedQuadrature(*TEST, *TRIAL );
    const TensorQuadrature* QUAD=dynamic_cast<const TensorQuadrature *>(quad.get());
    
    MacroPartition macros =makeMacros(TEST,TRIAL,QUAD,macroSize);
    MMatrix result=initMatrix(TEST,TRIAL,QUAD);
    IMatrix macroInf(domDim,macroInfS);
    
    for (int i=0;i<domDim;++i)
        {
        macroNum(i) = macros[i].size();
        macroCom(i) = i==0? 1 : macroCom(i-1)*macroNum(i-1);
        macroTST(i) = (*TEST)[i].size();
        macroTRL(i) = (*TRIAL)[i].size();
        }
    DoNotOptimize(result);
    auto t_end = std::chrono::high_resolution_clock::now();
    time_macro_setup.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    
    maxThread = std::min<unsigned int>(std::thread::hardware_concurrency(), maxThread);
    ThreadPool workers(maxThread);
    for (index_t r=0; r< (1<<domDim); ++r)
        {
        for (int c=0;c<domDim;++c)
            macroPos(c)= (r& (1<<c)) ? 1: 0;
        
        while ( macroPos.dot(macroCom) < macroNum.prod() )
            {
            workers.enqueue(assembleMacroelement, std::ref(eq), std::ref(macros), std::ref(result), macroInf );
            int c=0;
            for (;c<domDim;++c)
                {
                macroPos(c)+=2;
                if (macroPos(c)<macroNum(c) )
                    break;
                else
                    macroPos(c)= (r& (1<<c)) ? 1: 0;
                }
            if (c==domDim) break;
            }
                workers.waitAll();
        }
    return result;
    }
#undef macroInfS
#undef macroNum
#undef macroCom
#undef macroPos
#undef macroTST
#undef macroTRL
#undef macroSUB
}
