/** This file is part of the IgASF library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <map>
#include <tuple>
#include <json.hpp>
#include <iostream>
using Json = nlohmann::json;

template <typename Base>
struct Factory
{
    typedef Base* (*Maker) (const Json &j);
    struct SubType {
        std::string name;
        Maker       maker;
        } ;
    
    static  std::vector<SubType> m_subTypes;
    static  Base* makeSubType     (const Json &j)
    {
        std::string type=j["type"];
       auto st = std::find_if(m_subTypes.begin(),m_subTypes.end(), [type] (const SubType & a){return a.name==type;} );
        if (st!=m_subTypes.end())
            return st->maker(j);
        else
        {
            std::cerr<<"Subtype not known for \""<<type<<"\", returning nullptr"<<std::endl;
            std::cerr<<"Address is  "<<&m_subTypes<< std::endl;
            return nullptr;
        }

    }
    
    template <typename Derived>
    static  size_t registerSubType (std::string type, Maker  m = ConstructSubClass<Derived> )
    {
        auto st = std::find_if(m_subTypes.begin(),m_subTypes.end(), [type] (const SubType & a){return a.name==type;} );
//        std::cerr<<"Register "<<type<<std::endl;
//        std::cerr<<"Address is  "<<&m_subTypes<< std::endl;

        if (st!=m_subTypes.end())
        {
            std::cerr<<"Double subtype registration for \""<<type<<"\", request ignored"<<std::endl;
            std::cerr<<"Address is  "<<&m_subTypes<< std::endl;
            return -1;
        }
        else
        {
            m_subTypes.push_back({type,m});
            return m_subTypes.size();
        }
    }
    
    template <typename Derived>
    static  Base*  ConstructSubClass (const Json &j)
    {
        Derived * ret=new Derived();
        from_json(j,*ret);
        return ret;
    }

    virtual Json toJson() const =0;
};
template< typename Base > std::vector<typename Factory<Base>::SubType> Factory<Base>::m_subTypes;



template <typename Base, typename Derived, typename Middle=Base>
struct SubType : public Middle
{
    static size_t ID; // used by derived classes to hold the reference to their entry
    static Base * fromJson (const Json &j)
    {
        Derived * ret=new Derived();
        from_json(j,*ret);
        return ret;
    }
    virtual Json toJson() const
    {
        Json ret;
        to_json(ret, *dynamic_cast<const Derived*>(this));
        return ret;
    }
};

#define REGISTER_SUBTYPE( Base, Derived) \
template <>\
size_t SubType<Base, Derived>::ID=Base::registerSubType<Derived>( #Derived );

#define REGISTER_INDEX (Base) \
template <> \
struct Factory<Base>

// std::vector<Factory<Base>::SubType>  ::m_subTypes;


template <typename Base>
static void to_json(Json& j, const std::unique_ptr<Base>& p)
{ j=p->toJson(); }

template <typename Base>
void from_json(const Json& j, std::unique_ptr<Base>& p)
{ p=std::move(std::unique_ptr<Base>(Base::makeSubType(j))); }



