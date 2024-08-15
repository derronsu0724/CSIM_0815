/**
 * @file Ideal AC Voltage Source Model
 */

/*
 *  FastCSIM Copyright (C) 2021 cassuto                                    
 *  This project is free edition; you can redistribute it and/or           
 *  modify it under the terms of the GNU Lesser General Public             
 *  License(LGPL) as published by the Free Software Foundation; either      
 *  version 2.1 of the License, or (at your option) any later version.     
 *                                                                         
 *  This project is distributed in the hope that it will be useful,        
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU      
 *  Lesser General Public License for more details.                        
 */

#include <cmath>
#include "Constants.h"
#include "VAC.h"

namespace csimModel
{

    VAC::VAC(MODELBASE_CONSTRUCTOR_DEF)
        : ModelBase(MODELBASE_CONSTRUCTOR_VAR),
          m_Vp(0.0),
          m_omega(0.0),
          m_phase(0.0)
    {
        property().addProperty("Vp", Variant(Variant::VariantDouble).setDouble(5.0), nullptr, csimModel::PropertyBag::Required | csimModel::PropertyBag::Write | csimModel::PropertyBag::Read);
        property().addProperty("freq", Variant(Variant::VariantDouble).setDouble(50.0), nullptr, csimModel::PropertyBag::Required | csimModel::PropertyBag::Write | csimModel::PropertyBag::Read);
        property().addProperty("phase", Variant(Variant::VariantDouble).setDouble(0.0), nullptr, csimModel::PropertyBag::Write | csimModel::PropertyBag::Read);
    }

    VAC::~VAC()
    {
    }

    int VAC::configure()
    {
        m_Vp = property().getProperty("Vp").getDouble();
        m_omega = 2 * M_PI * property().getProperty("freq").getDouble();
        m_phase = M_PI * property().getProperty("phase").getDouble() / 180.0;

        m_E = MComplex(m_Vp * std::cos(m_phase), m_Vp * std::sin(m_phase));

        resizeModel(2, 0, 1);
        return 0;
    }

    int VAC::prepareDC()
    {
        return 0;
    }
    int VAC::prepareAC()
    {
        return 0;
    }
    int VAC::prepareTR()
    {
        return 0;
    }

    int VAC::iterateDC()
    {
        //std::cout << "getBranch(0):" <<getBranch(0)<< std::endl;
        //std::cout << "getNode(0):" <<getNode(0)<< std::endl;
        //std::cout << "getNode(0):" <<getNode(1)<< std::endl;
        unsigned int k = getBranch(0);
        addB(getNode(0), k, +1.0);
        addB(getNode(1), k, -1.0);
        addC(k, getNode(0), 1.0), addC(k, getNode(1), -1.0);
        addE(k, 0.0);
        return 0;
    }

    int VAC::iterateAC(double omega)
    {
        //std::cout << "getBranch(0):" <<getBranch(0)<< std::endl;
        //std::cout << "getNode(0):" <<getNode(0)<< std::endl;
        //std::cout << "getNode(0):" <<getNode(1)<< std::endl;
        unsigned int k = getBranch(0);
        addB(getNode(0), k, +1.0);
        addB(getNode(1), k, -1.0);
        addC(k, getNode(0), 1.0), addC(k, getNode(1), -1.0);
        addE(k, m_E);
        return 0;
    }

    int VAC::iterateTR(double tTime)
    {
        unsigned int k = getBranch(0);
        addB(getNode(0), k, +1.0);
        addB(getNode(1), k, -1.0);
        addC(k, getNode(0), 1.0), addC(k, getNode(1), -1.0);
        addE(k, m_Vp * std::sin(m_omega * tTime + m_phase));
        return 0;
    }

}

extern "C" csimModel::ModelBase *createModel(MODELBASE_CONSTRUCTOR_DEF)
{
    return new csimModel::VAC(MODELBASE_CONSTRUCTOR_VAR);
}

extern "C" void deleteModel(csimModel::ModelBase *model)
{
    delete model;
}

const PropertyMdlDescriptor *mdlDescriptors = nullptr;

const size_t numMdlDescriptors = 0;

const ModelDescriptor descriptor = {
    /* id */
    "VAC",
    /* description */
    "Generic ideal AC Voltage Source"};
