/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>

#include <maps/bc_map.h>
#include <maps/rational_map.h>

#include <bases/bspline.h>
#include <bases/tensor_basis.h>
#include <quadrature/tensor_quadrature.h>
#include <assembling/second_order.h>

struct iequal
{
    bool operator()(int c1, int c2) const
    {
        return std::toupper(c1) == std::toupper(c2);
    }
};

bool iequals(const std::string& str1, const std::string& str2)
{
    return str1.size() == str2.size()
        && std::equal(str1.begin(), str1.end(), str2.begin(), iequal());
}

inline std::string slurp (const std::string& path) {
    std::ostringstream buf; std::ifstream input (path.c_str()); buf << input.rdbuf(); return buf.str();
}

GeoPtr   getGeometry  (int dim ); // identity
BasisPtr getBasis     (std::vector<int> deg, std::vector<int> ele, std::vector<int> mul);
BasisPtr getBasis     (std::vector<int> deg, std::vector<int> ele, std::vector<int> mul, const GeoPtr &geo);

void print_help(char** args)
{
    std::cout<<"\n"<<args[0]
    <<" [-geo name] [-dim dimension] [-d degree] [-n element]\n"
    <<"           [-s smoothness] [-A coefs] [-b coefs] [-c coefs]\n\n"
    <<"Generates a test file (problem description file) and writes it to\n"
    <<"the standard output. In general, the parameters are optional and\n"
    <<"can be given in any order.\n\n"
    <<"  -geo name      name of domain geometry; see geometries.json for\n"
    <<"                 available geometries\n"
    <<"  -dim dimension specifies the domain to be hypercube of that dimension\n"
    <<"Either -geo or -dim must be given.\n\n"
    <<"  -d degree      polynomial degree of the test and trial functions; defaults to 1\n"
    <<"  -n element     number of polynomial segments for the test and trial spaces;\n"
    <<"                 defaults to 1\n"
    <<"  -s smoothness  smoothness of the test and trial spaces; defaults to degree-1\n"
    <<"Different values per direction are possible using -dd, -nn and -ss; this requires\n"
    <<"-geo or -dim to be specified before.\n\n"
    <<"The assembling procedure assembles for the problem\n"
    <<"  ∇u A ∇w + b.∇u w +c uw = fw\n"
    <<"The A, b and c options must follow -geo or -dim options and:\n"
    <<"  -A coefs       coefs must be dim^2 numbers, or ID for the identity;\n"
    <<"                 the numbers are in row major order: A(1,1) A(1,2) ... A(dim,dim);\n"
    <<"                 defaults to 0 ... 0.\n"
    <<"  -b coefs       coefs must be dim numbers; defaults to 0 ... 0\n"
    <<"  -c coef        coef is one number; defaults to 1\n"
    <<std::endl;
}

struct Options
{
    FMatrix             A; // Coefficients of PDE
    FMatrix             b; //  ---"---
    real_t              c; //  ---"---
    std::vector<int>    s; // smoothness
    std::vector<int>    d; // degree
    std::vector<int>    n; // number of elements
    int               dim;
    Json              geo; // geometry name
    
    Options() : c(1), dim(2) {}
};

int to_int(const char* in)
{
    if (!in)
        throw std::runtime_error("Reached end; expected integer.");
    try {
        return std::stoi(std::string(in));
        } 
    catch(...) {
        throw std::runtime_error(std::string("\"") + in + "\" cannot be interpreted as integer.");
        }
}

double to_double(const char* in)
{
    if (!in)
        throw std::runtime_error("Reached end; expected floating point number.");
    try {
        return std::stod(std::string(in));
        } 
    catch(...) {
        throw std::runtime_error(std::string("\"") + in + "\" cannot be interpreted as floating point number.");
        }
}

std::string to_string(const char* in)
{
    if (!in)
        throw std::runtime_error("Reached end; expected string.");
    return std::string(in);
}

Options parse_arguments(int argn, char** args)
{
    bool haveA=false;
    bool haveB=false;
    bool haveC=false;

    bool haveDim=false;
    bool haveEle=false;
    bool haveReg=false;
    bool haveDeg=false;
    bool haveSingleEle=true;
    bool haveSingleReg=true;
    bool haveSingleDeg=true;

    int n=1;
    int d=1;
    int s=INT_MAX; // very big so that it default to d-1

    Options result;
    char** end=args+argn;
    while (args<end)
        {
        if (iequals(std::string(args[0]),"-A"))
            {
            if ( haveA  ) throw std::runtime_error("Cannot give -A twice.");
            if (!haveDim) throw std::runtime_error("Domain must be known before -A.");
            ++args;
            result.A.resize(result.dim,result.dim);
            if (iequals(std::string(args[0]),"ID"))
                {
                ++args;
                result.A.setZero();
                for (int i=0;i<result.dim;++i)
                    result.A(i,i)=1;
                }
            else
                {
                for (int i=0;i<result.dim;++i)
                    for (int j=0;j<result.dim;++j)
                        result.A(i,j)=to_double(*(args++));
                }
            haveA=true;
            continue;
            }
        if(iequals(std::string(args[0]),"-b"))
            {
            if ( haveB  ) throw std::runtime_error("Cannot give -b twice.");
            if (!haveDim) throw std::runtime_error("Domain must be known before -b.");
            ++args;
            result.b.resize(result.dim,1);
            for (int j=0;j<result.dim;++j)
                result.b(j)=to_double(*(args++));
            haveB=true;
            continue;
            }
        if(iequals(std::string(args[0]),"-c"))
            {
            if ( haveC  ) throw std::runtime_error("Cannot give -c twice.");
            ++args;
            result.c=to_double(*(args++));
            haveC=true;
            continue;
            }
        if(iequals(std::string(args[0]),"-s"))
            {
            if ( haveReg) throw std::runtime_error("Cannot give -s or -ss twice.");
            ++args;
            s=to_int(*(args++));
            haveReg=true;
            continue;
            }
        if(iequals(std::string(args[0]),"-ss"))
            {
            if ( haveReg) throw std::runtime_error("Cannot give -s or -ss twice.");
            if (!haveDim) throw std::runtime_error("Domain must be known before -ss.");
            ++args;
            for (int i=0;i< result.dim;++i )
                result.s.push_back(to_int(*(args++)));
            haveReg=true;
            haveSingleReg=false;
            continue;
            }
        if(iequals(std::string(args[0]),"-d"))
            {
            if ( haveDeg) throw std::runtime_error("Cannot give -d or -dd twice.");
            ++args;
            d=to_int(*(args++));
            haveDeg=true;
            continue;
            }
        if(iequals(std::string(args[0]),"-dd"))
            {
            if ( haveDeg) throw std::runtime_error("Cannot give -d or -dd twice.");
            if (!haveDim) throw std::runtime_error("Domain must be known before -dd.");
            ++args;
            for (int i=0;i< result.dim;++i )
                result.d.push_back(to_int(*(args++)));
            haveDeg=true;
            haveSingleDeg=false;
            continue;
            }
        if(iequals(std::string(args[0]),"-n"))
            {
            if ( haveEle) throw std::runtime_error("Cannot give -n or -nn twice.");
            ++args;
            n=to_int(*(args++));
            haveEle=true;
            continue;
            }
        if(iequals(std::string(args[0]),"-nn"))
            {
            if ( haveEle) throw std::runtime_error("Cannot give -n or -nn twice.");
            if (!haveDim) throw std::runtime_error("Domain must be known before -nn.");
            ++args;
            for (int i=0;i< result.dim;++i )
                result.n.push_back(to_int(*(args++)));
            haveEle=true;
            haveSingleEle=false;
            continue;
            }
        if(iequals(std::string(args[0]),"-dim"))
            {
            if ( haveDim) throw std::runtime_error("Cannot give domain twice.");
            ++args;
            result.dim=to_int(*(args++));
            haveDim=true;
            continue;
            }
        if(iequals(std::string(args[0]),"-geo"))
            {
            if ( haveDim) throw std::runtime_error("Cannot give domain twice.");
            ++args;
            static std::string path("geometries.json");
            
            std::string text=slurp(path);
            Json data=Json::parse(text);
            std::string name=to_string(*(args++));
            try {
                GeoPtr geoMap=data[name].get<GeoPtr>();
                result.geo=geoMap;
                result.dim=geoMap->domainDim();
                haveDim=true;
                }
            catch (...) {
                throw std::runtime_error("Invalid geometry name.");
                }
            continue;
            }
        throw std::runtime_error(std::string("Unknown option \"") + args[0] + "\".");
        }
    if (!haveDim) throw std::runtime_error("Domain must be given.");
    if (haveSingleEle) result.n =std::vector<int>(result.dim,n);
    if (haveSingleDeg) result.d =std::vector<int>(result.dim,d);
    if (haveSingleReg) result.s =std::vector<int>(result.dim,s);

    return result;
}

int main (int argn, char** args)
{
    Options opt;
    try {
        opt = parse_arguments(argn-1, args+1);
        }
    catch(const std::exception& e) {
        std::cout << "\nThe following error occured: " << e.what() << std::endl;
        print_help(args);
        return 1;
        }

    GeoPtr geo;
    try {
        geo=opt.geo.get<GeoPtr>();
        }
    catch(...) {}

    if (geo.get()!=nullptr)
        opt.dim=geo->domainDim();
    else
        geo=getGeometry  (opt.dim);
    
    std::vector<int> mult(opt.d.size());
    for (size_t i=0; i< mult.size(); ++i)
        mult[i]=opt.d[i]-std::min(opt.s[i], opt.d[i]-1);
    Json result;
    
    BasisPtr bas=getBasis (opt.d,opt.n,mult,geo);
    QuadPtr  qua=getRecommendedQuadrature(static_cast<TensorBasis&>(*bas), static_cast<TensorBasis&>(*bas) );
    result["geometry"]   =geo;
    result["test"]       =bas;
    result["trial"]      =bas;
    result["quadrature"] =qua;
    result["EqCoefs"]    =EqCoef(opt.dim, opt.A,opt.b,opt.c);
    
    std::cout<<result<<std::endl;
    return 0;
}


FMatrix getGeoCoefs(int dim)
{
    int num=1<<dim; // = 2^dim
    FMatrix res(dim,num);
    for (int i=0; i< num; ++i)
        {
        for (int j=0;j<dim;++j)
            res(j,i)= i&(1<<j)? 1:0;
        }
    return res;
}

BasisPtr getBasis( std::vector<int> deg, std::vector<int> ele, std::vector<int> mul)
{
    int dim=deg.size();
    std::vector<BasisPtr> comps;
    for (int d=0;d<dim;++d)
    {
    std::vector<real_t> begs(deg[d]+1,0);
    std::vector<real_t> ends(deg[d]+1,ele[d]);
    std::vector<real_t> knts=begs;

    for (int i=1;i<ele[d];++i)
        for (int j=0;j<mul[d];++j)
            knts.push_back(i);
    knts.insert(knts.end(),ends.begin(),ends.end());
    comps.push_back(BasisPtr(new Bspline(deg[d],knts)));
    }
    return BasisPtr(new TensorBasis(comps));
}

BasisPtr getBasis(std::vector<int> deg, std::vector<int> ele, std::vector<int> mul, const GeoPtr &geo)
{
    int dim=geo->domainDim();
    
    const TensorBasis *bas;
    const BasisCoefficientMap* map;
    const RationalMap* rmap=dynamic_cast<const RationalMap*>(geo.get());
    if (rmap)
        map=dynamic_cast<const BasisCoefficientMap*>(rmap->underlying().get());
    else
        map=dynamic_cast<const BasisCoefficientMap*>(geo.get());
    
    if (map) bas=dynamic_cast<const TensorBasis*>(map->m_basis.get());
    else throw std::runtime_error("only works with NURBS or Bspline maps");
    
    std::vector<BasisPtr> comps;
    for (int d=0;d<dim;++d)
    {
        const std::vector<real_t> &g_knt=dynamic_cast<const Bspline&>((*bas)[d]).knots();
        const int addMult=std::max(0, deg[d]-dynamic_cast<const Bspline&>((*bas)[d]).degree());

        real_t st=g_knt.front();
        real_t en=g_knt.back();
        real_t step=(en-st)/real_t(ele[d]);
    
        std::vector<real_t> begs(deg[d]+1,st);
        std::vector<real_t> ends(deg[d]+1,en);
        std::vector<real_t> knts=begs;
        std::vector<real_t> merged;

        for (int i=1;i<ele[d];++i)
        for (int j=0;j<mul[d];++j)
        {
            knts.push_back(i*step+st);
        }
        knts.insert(knts.end(),ends.begin(),ends.end());
        auto curK=knts.begin();
        auto endK=knts.end();
        auto curG=g_knt.begin();
        auto endG=g_knt.end();
        while(curG<endG && curK<endK)
        {
            if (*curK >= *curG)
            {
                auto nextG=std::find_if(curG,endG,[curG](auto v){return v!=*curG;});
                auto nextK=std::find_if(curK,endK,[curG](auto v){return v!=*curG;});
                int  mult=std::max(nextG-curG+addMult, nextK-curK);
                for (int i=0; i<mult;++i)
                    merged.push_back(*curG);
                curG=nextG;
                curK=nextK;
                continue;
            }
            if (*curK < *curG) {merged.push_back(*curK); ++curK; }
        }
        merged.insert(merged.end(),curK,endK);
        merged.insert(merged.end(),curG,endG);
        comps.push_back(BasisPtr(new Bspline(deg[d],merged)));
    }
    return BasisPtr(new TensorBasis(comps));
}

GeoPtr getGeometry(int dim)
{
    std::vector<int> data(dim,1);

    BasisPtr bas  = getBasis(data, data, data);
    FMatrix  coef = getGeoCoefs(dim);
    return GeoPtr(new BasisCoefficientMap(bas, coef));
}



