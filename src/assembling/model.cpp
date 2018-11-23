/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <map>
#include <set>
#include <iostream>

#include <assembling/model.h>
#include <assembling/recursive_nD.h>
#include <assembling/apply.h>
#include <maps/transform_coefs.h>
#include <tools/timing.h>

void Model::assemble(const TensorBasis &test, const TensorBasis &trial, const TensorQuadrature  &quad, MMatrix &output) const
{
    int dim = quad.domainDim();
    
    // coefficients
    auto t_start = std::chrono::high_resolution_clock::now();
    std::vector<Part> data=this->initParts(quad);
    auto t_end = std::chrono::high_resolution_clock::now();
    time_eval_coef.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    
    // basis functions
    t_start = std::chrono::high_resolution_clock::now();
    auto tsts=test.evaluateComponents (getTestRequest(dim,data),quad.grid());
    auto trls=trial.evaluateComponents(getTrialRequest(dim,data),quad.grid());
    quad.applyToValues(tsts);
    t_end = std::chrono::high_resolution_clock::now();
    time_eval_bases.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    
    // elements
    t_start = std::chrono::high_resolution_clock::now();
    auto elesData=computeElementSplitting(test,trial,quad);
    std::vector<view<const index_t>>  eles(dim);
    for (int i=0;i<dim;++i)eles[i]=view<const index_t>(elesData[i]);
    
    // sparsity
    std::vector<Sparsity>    sprs = getBilSparsities(tsts, trls);
    std::vector<Sparsity>    krns = getKroSparsities(sprs);
    
    // temporaries
    index_t tmp_size=0;
    for (int i=0; i<dim-1;++i)
        tmp_size+=krns[i].nnzs();
    static thread_local std::vector<real_t>          memsD;
    static thread_local std::vector<view<real_t>>    memsV;
    static thread_local std::vector<rview<2,real_t>> memsW;
    memsD.resize(MAX_TMP*tmp_size);
    memsV.resize(MAX_TMP*(dim-1));
    memsW.resize(dim-1);
    
    index_t start=0;
    for (int c=0; c<dim-1;++c)
        {
        index_t sz=krns[c].nnzs();
        for (int r=0; r<MAX_TMP;++r)
            {
            memsV[c*MAX_TMP+r]=view<real_t>(memsD).subView(start,sz);
            start+=sz;
            }
        memsW[c]=rview<2,real_t>(memsV).subView(c*MAX_TMP,MAX_TMP);
        }
    
    // init output
    output.coefs.resize(krns.back().nnzs());
    view<real_t>(output.coefs).vector().setZero();
    t_end = std::chrono::high_resolution_clock::now();
    time_compute_structure.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());

    
    // Computation
    t_start = std::chrono::high_resolution_clock::now();
    for (auto &entry: data)
        {
        recursiveAssemble(
                          view<const BasisValues>(tsts),
                          view<const BasisValues>(trls),
                          entry.test,
                          entry.trial,
                          entry.coefs,
                          view<real_t> (output.coefs),
                          view<const Sparsity>(sprs),
                          view<const Sparsity>(krns),
                          rview<2,const index_t>(eles),
                          rview<3,real_t>(memsW)
                          );
        }
    std::swap<Sparsity>(krns.back(), output);
    t_end = std::chrono::high_resolution_clock::now();
    time_assemble.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    this->done();
}

void Model::apply   (const TensorBasis &test, const TensorBasis &trial, const TensorQuadrature  &quad, view<const real_t> in_v, view<real_t> out_v) const
{
    int dim = quad.domainDim();
    
    // coefficients
    auto t_start = std::chrono::high_resolution_clock::now();
    std::vector<Part> data=this->initParts(quad);
    auto t_end = std::chrono::high_resolution_clock::now();
    time_eval_coef.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    
    // basis functions
    t_start = std::chrono::high_resolution_clock::now();
    auto tsts=test.evaluateComponents (getTestRequest(dim,data),quad.grid());
    auto trls=trial.evaluateComponents(getTrialRequest(dim,data),quad.grid());
    quad.applyToValues(tsts);
    t_end = std::chrono::high_resolution_clock::now();
    time_eval_bases.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    
    // elements
    t_start = std::chrono::high_resolution_clock::now();
    auto elesData=computeElementSplitting(test,trial,quad);
    std::vector<view<const index_t>>  eles(dim);
    for (int i=0;i<dim;++i)eles[i]=view<const index_t>(elesData[i]);
    
    // init output and tmps
    index_t max_size=0;
    {
    index_t size=in_v.size();
    for (int i=dim-1; i>=0;--i)
        {
        size/=trls[i].rows();
        size*=trls[i].cols();
        max_size=std::max(max_size, size);
        }
    }
    {
    index_t size=quad.size();
    for (int i=dim-1; i>=0;--i)
        {
        size/=tsts[i].rows();
        size*=tsts[i].cols();
        max_size=std::max(max_size, size);
        }
    }
    
    OwningView<real_t> tmp(max_size);
    OwningView<real_t> eval_mem(quad.size());
    OwningView<real_t> int_mem(out_v.size());

    out_v.vector().setZero();
    t_end = std::chrono::high_resolution_clock::now();
    time_compute_structure.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());

    
    // Computation
    t_start = std::chrono::high_resolution_clock::now();
    for (auto &entry: data)
        {
        KroneckerApply(
                          view<const BasisValues>(tsts),
                          view<const BasisValues>(trls),
                          entry.test,
                          entry.trial,
                          entry.coefs,
                          eval_mem,
                          int_mem,
                          tmp,
                          in_v,
                          out_v
                          );
        }
    t_end = std::chrono::high_resolution_clock::now();
    time_assemble.fetch_add(
        std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count());
    this->done();
}

std::vector<Sparsity> Model::getBilSparsities(const std::vector<BasisValues> &tsts, const std::vector<BasisValues> &trls)
{
    int dim=tsts.size();
    std::vector<Sparsity>    sprs(dim);
    for (int c=0; c<dim;++c)
        sprs[c]=bilinearSparsity(tsts[c],trls[c]);
    return sprs;
}

std::vector<Sparsity> Model::getKroSparsities(const std::vector<Sparsity> &sprs)
{
    int dim=sprs.size();
    std::vector<Sparsity>    krns(dim);
    krns[0]=Sparsity(sprs[0]);
    for (int c=1; c<dim;++c)
        krns[c]=kroneckerSparsity(sprs[c],krns[c-1]);
    return krns;
}

std::vector<std::vector<index_t>> Model::computeElementSplitting(const TensorBasis &test, const TensorBasis &trial, const TensorQuadrature  &quad)const
{
    int dim=test.domainDim();
    std::vector<std::vector<index_t>> elesData(dim);
    for (int i=0;i<dim;++i)
        {
        auto &quadC=quad[i];
        auto &nodes=quadC->nodes();
        
        auto &testBrk = test[i].breaks();
        auto begTest= ++(testBrk.begin());
        auto endTest= testBrk.end();
        auto &trialBrk = trial[i].breaks();
        auto begTrial = ++( trialBrk.begin());
        auto endTrial = trialBrk.end();
        
        index_t pt=0;
        while (begTest !=endTest && begTrial !=endTrial)
            {
            while(pt< quadC->size() && nodes(0,pt) < *begTest &&  nodes(0,pt) < *begTrial)
                ++pt;
            elesData[i].push_back(pt);
            if (pt==quadC->size() ) break;
            while ( begTest  !=endTest  && nodes(0,pt) >= *begTest )  ++begTest;
            while ( begTrial !=endTrial && nodes(0,pt) >= *begTrial ) ++begTrial;
            }
        }
    return elesData;
}

std::vector<PartialDerivative> Model::getModelRequest(Role role, int component,const std::vector<Part> &data )const
{
    std::set<PartialDerivative> res;
    for (auto &entry:data)
        { res.insert((entry.*role)[component]); }
    return std::vector<PartialDerivative>(res.begin(), res.end());
}

TensorBasis::ComponentRequest Model::getTestRequest (int dim,const std::vector<Part> &data)const
{
    TensorBasis::ComponentRequest res(dim);
    for (int i=0; i<dim;++i)
        res[i]=getModelRequest(&Part::test, i,data);
    return res;
}

TensorBasis::ComponentRequest Model::getTrialRequest( int dim,const std::vector<Part> &data)const
{
    TensorBasis::ComponentRequest res(dim);
    for (int i=0; i<dim;++i)
        res[i]=getModelRequest(&Part::trial, i,data);
    return res;
}

void to_json  (Json& j, const Model::Part& p)
{
    j["tstD"]=p.test;
    j["trlD"]=p.trial;
    j["vals"]=p.coefs;
}
