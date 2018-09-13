/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <iostream>
#include <fstream>
#include <chrono>
#include <exception>
#include <algorithm>
#include <string>

#include <assembling/second_order.h>
#include <assembling/macroelement.h>

#include <tools/timing.h>
#include <tools/matrixio.h>
#include <json.hpp>

using Json = nlohmann::json;


struct problemData
{
    std::string       problem_file;
    Json              problem;
    std::string       m;
    bool              global;
    std::vector<int>  size;
    unsigned int      threads;
    std::string       output;
    std::string       log;
    problemData() : m("global"), global(true), threads(1), output(), log() {}
};

void print_help(char** args)
{
    std::cout<<"\n"<<args[0]
    <<" test_file [-o file] [-l file] [-m method] [-threads num]\n\n"
    <<"  test_file      the test file (as created with generateTest); use \"stdin\"\n"
    <<"                 to read from console.\n"
    <<"  -o file        specifies a file to write the matrix into; use \"stdout\"\n"
    <<"                 to write to console.\n"
    <<"  -l file        specifies a file to write log data into (append to).\n"
    <<"  -m method      method can be\n"
    <<"                    global     for global sum factorization (default)\n"
    <<"                    element    for element-wise sum factorization;\n"
    <<"                               this is the same as -m macro 1 ... 1 \n"
    <<"                    macroS     for macro elements of size of degree;\n"
    <<"                               this is the same as -m macro p[1] ... p[d]\n"
    <<"                    macroN     this is the same as -m macro p[1] ... p[d-1] 1\n"
    <<"                    macroR     this is the same as -m macro 1 p[2] ... p[d]\n"
    <<"                    macro s    this is the same as -m macro s ... s\n"
    <<"                    macro s1 ... sd\n" 
    <<"                               for macro elements of specified size; a\n"
    <<"                               value of -1 is replaced by the spline degree\n"
    <<"                               of the respective direction.\n"
    <<"  -threads num   to use num (default: 1) threads; has no effect for -m global.\n"
    <<std::endl;
}

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

std::string to_string(const char* in)
{
    if (!in)
        throw std::runtime_error("Reached end; expected string.");
    return std::string(in);
}

bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

problemData parse_args(int argn, char** args)
{
    char** end=args+argn;

    bool haveO = false;
    bool haveL = false;
    bool haveM = false;
    bool haveT = false;

    problemData r;

    if (argn<1) throw std::runtime_error("No test file is given.");

    // input
    r.problem_file = args[0];
    if (std::string(args[0])=="stdin")
        std::cin >> r.problem;
    else
        std::ifstream(args[0])>>r.problem;

    args++;
    
    while (args<end)
        {
        if (iequals(args[0],"-o"))
            {
            if ( haveO  ) throw std::runtime_error("Cannot give -o twice.");
            ++args;
            std::string fn = to_string(*(args++));
            r.output = fn;
            haveO = true;
            continue;
            }
        if (iequals(args[0],"-l"))
            {
            if ( haveL  ) throw std::runtime_error("Cannot give -l twice.");
            ++args;
            r.log = to_string(*(args++));
            haveL = true;
            continue;
            }
        if (iequals(args[0],"-m"))
            {
            if ( haveM  ) throw std::runtime_error("Cannot give -m twice.");
            ++args;
            r.m = to_string(*(args++));
            if (r.m=="global")
                r.global = true;
            else if (r.m=="element")
                r.global = false;
            else if (r.m=="macroS")
                r.global = false;
            else if (r.m=="macroN")
                r.global = false;
            else if (r.m=="macroR")
                r.global = false;
            else if (r.m=="macro")
            {
                r.global = false;
                try {
                    while (true)
                        {
                        const int sz = to_int(*(args++));
                        r.size.push_back(sz);
                        r.m += ":";
                        r.m += std::to_string(sz);
                        }
                    }
                catch (...) {}
                --args;
                if (r.size.size()==0)
                    throw std::runtime_error("-m macro requires valid sizes to be given.");
            }
            else
                throw std::runtime_error(std::string("Unknown method \"") + r.m + "\".");
            haveM = true;
            continue;
            }
        if (iequals(args[0],"-threads"))
            {
            if ( haveT  ) throw std::runtime_error("Cannot give -threads twice.");
            ++args;
            r.threads = to_int(*(args++));
            if (r.threads<1) throw std::runtime_error("Need at least 1 thread.");
            haveT = true;
            continue;
            }
        throw std::runtime_error(std::string("Unknown option \"") + args[0] + "\".");
        }
    if (r.threads>1 &&  r.m != "macroS")
        throw std::runtime_error(std::string("Parallel implementation is only available for \"macroS\"."));
    return r;
}

int main (int argn, char ** args)
{
    problemData data;
    try {
        data=parse_args(argn-1,args+1);
        }
    catch(std::exception& e) {
        std::cout << "\nThe following error occured: " << e.what() << std::endl;
        print_help(args);
        return 1;
    }
        
    // convert json to internal representation
    EqCoef eq;
    BasisPtr tstSp;
    BasisPtr trlSp;
    GeoPtr   geo;
    QuadPtr  quad;
    
    try {
        try { // eq is optional
            eq = data.problem.at("EqCoefs").get<EqCoef>();
            }
        catch(...){}
        tstSp = data.problem.at("test").get<BasisPtr>();
        trlSp = data.problem.at("trial").get<BasisPtr>();
        geo   = data.problem.at("geometry").get<GeoPtr>();
        quad  = data.problem.at("quadrature").get<QuadPtr>();
        }
    catch (const std::string &s) {
        std::cerr<<"Init abort with error: "<<s<<std::endl;
        return 1;
        }
    
    float  real_time, cpus_time;
    MMatrix res;
    
    try {
        auto real_start = std::chrono::high_resolution_clock::now();
        auto cpus_start = std::clock();
        SecondOrderModel m(eq, geo.get());
        TensorBasis &tstB=*dynamic_cast<TensorBasis*>(tstSp.get());
        TensorBasis &trlB=*dynamic_cast<TensorBasis*>(trlSp.get());
        QuadPtr quad=getRecommendedQuadrature(tstB,trlB);

        if (data.global)
            m.assemble(tstB,trlB,*static_cast<TensorQuadrature*>(quad.get()), res);
        else
            {
            if (data.m == "element")
                {
                data.size.resize(tstSp->domainDim());
                view<index_t>(data.size).vector().setOnes();
                }
            else if (data.m == "macroS")
                {
                data.size.resize(tstSp->domainDim());
                view<index_t>(data.size).vector().setConstant(-1);
                }
            else if (data.m == "macroN")
                {
                data.size.resize(tstSp->domainDim());
                view<index_t>(data.size).vector().setConstant(-1);
                data.size.back() = 1;
                }
            else if (data.m == "macroR")
                {
                data.size.resize(tstSp->domainDim());
                view<index_t>(data.size).vector().setConstant(-1);
                data.size[0] = 1;
                }
            res = macroelement::assembleParallel (m, tstSp, trlSp, data.threads, view<index_t>(data.size).vector());
            }

        auto real_end = std::chrono::high_resolution_clock::now();
        cpus_time = float(std::clock()-cpus_start) / CLOCKS_PER_SEC;
        real_time = std::chrono::duration<float>(real_end-real_start).count();
        }
    catch (std::exception& e) {
        std::cout << "An error occured during assembling: " << e.what() << "\n" << std::endl;
        return 1;
    }

    if (data.output!="stdout")
        std::cout << "Have assembled the corresponding " << res.rows() << "x" << res.cols() << " Galerkin matrix using "
            <<"sum factorization ("<<data.m<<").\n"
            <<"Time:\n"
            <<"    kron-struct: "<<float(time_compute_structure)/1000000<<"s\n"
            <<"    bases-eval:  "<<float(time_eval_bases       )/1000000<<"s\n"
            <<"    coefs-eval:  "<<float(time_eval_coef        )/1000000<<"s\n"
            <<"     -geo-eval:  "<<float(time_geo_compute      )/1000000<<"s\n"
            <<"     -geo-tran:  "<<float(time_geo_transform    )/1000000<<"s\n"
            <<"    assemble:    "<<float(time_assemble         )/1000000<<"s\n"
            <<"    macro-setup: "<<float(time_macro_setup      )/1000000<<"s\n"
            <<"    macro-add:   "<<float(time_add_macro        )/1000000<<"s\n\n"
            <<"Total time:      "<<real_time<<"s (real)\n"
            <<"                 "<<cpus_time<<"s (cpus)\n"
            <<std::endl;
    
    if (!data.output.empty())
        {
        if (data.output!="stdout")
            std::cout << "The assembled matrix is written to " << data.output << ".\n" << std::endl;
        outData output( const_cast<const MMatrix&>(res).matrix() );
        writeMatrix(output,data.output);
        }

    if (!data.log.empty())
        {
        if (data.output!="stdout")
            std::cout << "Log data is written to " << data.log << ".\n" << std::endl;

        bool exists = file_exists(data.log);

        std::fstream fs;
        fs.open(data.log, std::fstream::in | std::fstream::out | std::fstream::app);
        if (!fs.is_open())
            std::cout << "Failed to open file.\n" << std::endl;
        else
            {
            if (!exists)
                fs << "TestName\tMethod\tTotalClockTime\tTotalCpuTime\tPartialSparsity\tPartialBases\t"
                   << "PartialCoefficients\tPartialGeometry\tPartialTransformation\tPartialSumFactorization\t"
                   << "PartialMacroElementSetup\tPartialMacroElementMerge\n";

            fs << data.problem_file << '\t'
                << data.m << '\t'
                << real_time << '\t'
                << cpus_time << '\t'
                << float(time_compute_structure)/1000000 << '\t'
                << float(time_eval_bases       )/1000000 << '\t'
                << float(time_eval_coef        )/1000000 << '\t'
                << float(time_geo_compute      )/1000000 << '\t'
                << float(time_geo_transform    )/1000000 << '\t'
                << float(time_assemble         )/1000000 << '\t'
                << float(time_macro_setup      )/1000000 << '\t'
                << float(time_add_macro        )/1000000 << '\n' << std::flush;
            fs.close();
            }
        }

    return 0;
}


