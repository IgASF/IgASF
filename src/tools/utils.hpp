/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

template <typename T, typename U=T> // returns n(n-1)...(s+1)s
constexpr T factorial(T n, U s=U(1)) {
  return n>s ? (n * factorial(n - 1,s)) : 1;
}

template <typename T, typename U> // returns n(n-1)...(d+1) / [d(d-1)...2]
constexpr T binomial(T n, U d) {
  return factorial(n,d)/factorial(d);
}



template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T> &v)
{
  out<< "Vector: "<<v.size()
    << (v.size()? "(\n" : "(" );
  size_t c=0;
  for (auto t : v)
    out<<"["<<c++<<"]: "<<t<<"\n";
  out<<")\n\n";
  return out;
}


 template< typename FwdIter, typename Func >
 Func for_each_pair( FwdIter iterStart, FwdIter iterEnd, Func func )
 {
     if( iterStart == iterEnd )
        return func;

     FwdIter iterNext = iterStart;
     ++iterNext;

     for( ; iterNext != iterEnd; ++iterStart, ++iterNext )
     {
          func( *iterStart, *iterNext );
     }
     return func;
}

template <typename T, typename ... Ts>
auto min(T a, Ts... args)
{
    if constexpr (sizeof...(Ts)==0)
        return a;
    else
        return a<= min(args...) ? a : min(args...) ;
}

template <unsigned int N, template <int> typename F, typename ...ARGS>
auto dispatch( unsigned int n, ARGS... arg)
    {
        if (n==N) return F<N>()(arg...);
        if constexpr (N>1)
        {
            return dispatch<N-1, F>(n, arg...);
        }
    }

template <unsigned int N1, template <int, int> typename F>
struct Dispatch2Wrapper
{
    template <int M>
    using FF=F<N1,M>;
};

template <unsigned int N1,unsigned  int N2, template <int, int> typename F, typename ...ARGS>
auto dispatch2( unsigned int n1, unsigned int n2, ARGS... arg)
    {

    
        if (n1==N1) return dispatch<N2,Dispatch2Wrapper<N1,F>::template FF>(n2,arg...);
        if constexpr (N1>1)
        {
            return dispatch2<N1-1,N2,F>(n1,n2, arg...);
        }
    }

