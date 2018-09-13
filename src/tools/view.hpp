/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <algebra/eigen.h>

template <typename T>
struct view
{
protected:
    T*   m_beg;
    T*   m_end;
private:
    typedef typename std::remove_const<T>::type plainT;
    template <int R,int C, int O>
    using NonConstMatrix = Eigen::Map<Eigen::Matrix<plainT,R,C, O>>;
    template <int R,int C, int O>
    using ConstMatrix =  Eigen::Map<const Eigen::Matrix<plainT,R,C,O>>;
public:
    template <int R=-1,int C=-1,int O=0>
    using Matrix = typename std::conditional<std::is_const<T>::value, ConstMatrix<R,C,O>, NonConstMatrix<R,C,O>>::type;
    typedef Matrix<-1,1> Vector;
    template <int R=-1,int C=-1,int O=0>
    using Array = decltype(std::declval<Matrix<R,C,O>>().array());
public:
    view()                               :  m_beg(nullptr), m_end(nullptr){}
    view(T* beg, size_t size)            :  m_beg(beg),        m_end(beg+size) {}
    view(      std::vector<plainT> &vec) :  m_beg(vec.data()), m_end(vec.data()+vec.size()) {}
    view(const std::vector<plainT> &vec) :  m_beg(vec.data()), m_end(vec.data()+vec.size()) {}

    template <int i,int j, int o>
    view(Eigen::Matrix<plainT,i,j,o> &mat) :  m_beg(mat.data()), m_end(mat.data()+mat.size()) {}
    
    template <int i,int j, int o>
    view(const Eigen::Matrix<plainT,i,j,o> &mat) :  m_beg(mat.data()), m_end(mat.data()+mat.size()) {}


// access
    size_t   size()  const {return m_end-m_beg;}
    T*       begin() const {return m_beg;}
    T*       end()   const {return m_end;}
    T& operator[](index_t p) const {return m_beg[p];}

    T& back()  const   {return *(m_end-1);}
    T& front() const   {return *m_beg;}

    view subView(index_t start, index_t size) const { return view(m_beg+start, size); }
    view dropBack(index_t n=1)  const {return view(m_beg, size()-n);}
    view dropFront(index_t n=1) const {return view(m_beg+n, size()-n);}

    void popFront(index_t n=1) {m_beg+=n;}
    void popBack (index_t n=1) {m_end-=n;}
   
    // conversion
    operator view<const plainT> ()const {return view<const plainT>(m_beg, size());}
    // view as algebraic objects
    static constexpr int DO(int R, int ){return R==1? Eigen::RowMajor : 0;}
    
    template <int R,int C,int O=DO(R,C)>
    Matrix<R,C,O>  matrix(index_t r=R,index_t c=C)                     const {return Matrix<R,C,O>(m_beg,r,c); }
    Matrix<>     matrix(index_t r, index_t c) const {return Matrix<>(m_beg,r,c);    }
    Vector       vector()                     const {return Vector(m_beg,size());   }
    template <int R,int C,int O=DO(R,C)>
    Array<R,C,O>   array(index_t r=R,index_t c=C)                      const {return matrix<R,C,O>(r,c).array();  }
    Array<>      array(index_t r, index_t c)  const {return matrix(r,c).array();    }
    Array<-1,1>  varray()                     const {return vector().array();       }
    // change algebraic objects in place (by changing the view of the data, but leaving the data intact)
    template <int R,int C,int O=DO(R,C)>
    void  matrix(Matrix<R,C,O> &dst)                   const {new (&dst) Matrix<R,C,O>(m_beg,R,C); }
    void  matrix(index_t r, index_t c,Matrix<> &dst) const {new (&dst) Matrix<>(m_beg,r,c);    }
};

template <typename T>
struct OwningView : public view<T>
{
    using view<T>::m_beg;
    using view<T>::m_end;
    using view<T>::begin;
    using view<T>::end;

    OwningView(): view<T>(nullptr,0) {}
    OwningView( size_t size)  {  m_beg= size>0 ? new T[size]:nullptr; m_end= size>0 ?m_beg+size:nullptr;}

    OwningView( OwningView && other) { swap(std::move(other));}
    OwningView( const OwningView & other) : view<T>(new T[other.size()], other.size()) {std::copy(other.begin(),other.end(),begin());}
    template <typename It>
    OwningView( It beg, It end) : view<T>(new T[end-beg], end-beg) {std::copy(beg,end, begin());}
    
    OwningView& operator=(const OwningView&)=delete ;
    OwningView& operator=(OwningView&& other) { swap(other); return *this;}

    void swap(OwningView && other) {T *t1=m_beg,*t2=m_end; m_beg=other.m_beg; m_end=other.m_end; other.m_beg=t1; other.m_end=t2;}
    void swap(OwningView & other) {T *t1=m_beg,*t2=m_end; m_beg=other.m_beg; m_end=other.m_end; other.m_beg=t1; other.m_end=t2;}
    ~OwningView() {delete[] view<T>::m_beg;}
};

template <typename T>
void to_json  (Json& j, const view<T>& p)
{
    for (size_t i=0;i<p.size(); ++i)
    j[i]=p[i];
}

template <int n, typename R>
struct rview_impl
{
    typedef view<typename rview_impl<n-1,R>::type> type;
};

template <typename R>
struct rview_impl<0,R>
{
    typedef R type;
};

template <int n, typename R>
using rview=typename rview_impl<n,R>::type;
