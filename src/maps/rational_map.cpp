/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <maps/rational_map.h>
#include <tools/view.hpp>

namespace{
// implementation forward declarations
FMatrix evalRational  ( FMatrix &&tmp);
FMatrix evalRational  ( FMatrix &&val, FMatrix &&jac);
} // anonymous namespace


// public member functions

RationalMap::RationalMap(GeoPtr & geo)
: m_original(std::move(geo))
{}

const GeoPtr&  RationalMap::underlying() const
{return m_original;}


int RationalMap::targetDim() const
{return m_original->targetDim()-1;}

int RationalMap::domainDim() const
{return m_original->domainDim();}

FMatrix RationalMap::evaluate(const view<real_t> &points) const
{ return evalRational(m_original->evaluate(points)); }

FMatrix RationalMap::evaluate(const CartesianGrid &points) const
{ return evalRational(m_original->evaluate(points)); }


FMatrix RationalMap::jacobian(const view<real_t> &points) const
{ return evalRational(m_original->evaluate(points),m_original->jacobian(points)); }

FMatrix RationalMap::jacobian(const CartesianGrid &points) const
{ return evalRational(m_original->evaluate(points),m_original->jacobian(points)); }


// public IO functions
void to_json(Json& j, const RationalMap& p)
{ j = Json{ { "type", "RationalMap"}, {"original", p.underlying() } }; }

void from_json(const Json& j, RationalMap& p)
{
    if (j["type"]!="RationalMap") throw std::runtime_error("expected RationalMap, got "+j.at("type").get<std::string>() );
    if (p.underlying()!=nullptr)
        p.~RationalMap();
    GeoPtr ptr;
    from_json(j["original"],ptr);
    new (&p) RationalMap(ptr); // j["original"].get<GeoPtr>()
}


namespace {
// implementation

FMatrix evalRational (FMatrix &&tmp)
{ return tmp.bottomRows(1).asDiagonal().inverse()*tmp.topRows(tmp.rows()-1); }


template <int dom, int tar>
void pointLoop(index_t numPts, view<const real_t> val, view<const real_t> jac, view<real_t> out)
{
    for (index_t p=0; p<numPts; ++p)
        {
        auto J=jac.subView(p*dom*(tar+1),dom*(tar+1)).matrix<tar+1,dom>();
        auto V=val.subView(p*(tar+1),tar).matrix<tar,1>();
        const real_t den=val[(tar+1)*p+tar];
        
        out.subView(p*dom*tar,dom*tar).matrix<tar,dom>()
        = ( J.template topRows<tar>()*den
           -V*J.template bottomRows<1>())/den/den;
        //        if (std::getenv("debug")) { std::cout << V << "\n\n"<<J <<"\n\n"<< den<<"\n\n"<<out.subView(p*dom*tar,dom*tar).matrix<tar,dom>()<<"\n\n"<<std::endl; }
        
        }
}

FMatrix evalRational (FMatrix &&val, FMatrix &&jac)
{
    //    if (std::getenv("debug")) { std::cout << val << "\n\n "<<jac <<"\n"<<std::endl; }
    
    int tar=val.rows()-1;
    int dom=jac.rows()/(tar+1);
    int numPts=val.cols();
    FMatrix result(dom*tar,numPts);
    view<real_t> res(result);
    view<const real_t> valV(val);
    view<const real_t> jacV(jac);
    switch ( (dom-1)*4+(tar-1) )
    {
        // dom=1
        case 0: pointLoop<1,1>(numPts,valV,jacV,res); break;
        case 1: pointLoop<1,2>(numPts,valV,jacV,res); break;
        case 2: pointLoop<1,3>(numPts,valV,jacV,res); break;
        case 3: pointLoop<1,4>(numPts,valV,jacV,res); break;
        // dom=2
        case 5: pointLoop<2,2>(numPts,valV,jacV,res); break;
        case 6: pointLoop<2,3>(numPts,valV,jacV,res); break;
        case 7: pointLoop<2,4>(numPts,valV,jacV,res); break;
        // dom=3
        case 10: pointLoop<3,3>(numPts,valV,jacV,res); break;
        case 11: pointLoop<3,4>(numPts,valV,jacV,res); break;
        // dom 4
        case 15: pointLoop<4,4>(numPts,valV,jacV,res); break;
        default:
        assert(false);
    }
    //    if (std::getenv("debug")) { std::cout << result << "\n\n " <<"\n"<<std::endl; }
    return result;
}

} // anonymous namespace

