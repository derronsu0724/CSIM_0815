/**
 * @file DC analysis
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

#include "csim/model/ModelBase.h"
#include "csim/utils/errors.h"
#include "csim/internal/Netlist.h"
#include "csim/internal/Circuit.h"
#include "csim/internal/Dataset.h"
#include "AnalyzerOP.h"
#include "log.h"
namespace csim
{

    AnalyzerOP::AnalyzerOP(Circuit *circuit)
        : AnalyzerBase(circuit)
    {
    }
    AnalyzerOP::~AnalyzerOP()
    {
    }

    int AnalyzerOP::analyze(Dataset *dataset)
    {
        /*
         * Format data set
         */
        dataset->clear();
        dataset->setName("DC quiescent operating point analysis");
        unsigned int N = getNumInterestNodes(), M = getNumInterestBranches();
        std::vector<Variable *> dvolt(N), dcurrent(M);
        for (unsigned int i = 0; i < N; ++i)
        {
            std::string varName = makeVarName("V", getInterestNode(i));
            dvolt[i] = &dataset->addDependentVar("voltage", varName);
        }
        for (unsigned int i = 0; i < M; ++i)
        {
            std::string varName = makeVarName("I", getInterestBranch(i));
            dcurrent[i] = &dataset->addDependentVar("current", varName);
        }

        UPDATE_RC(circuit()->initMNA());
#if defined(ENABLE_SPICE_COMPATIBLE)
        UPDATE_RC(circuit()->spiceCompatible()->setFlagsOP());
#endif
        UPDATE_RC(circuit()->prepareMNA(this));
        UPDATE_RC(circuit()->solveMNA(this));
        UPDATE_RC(circuit()->saveOP());

        /* Save the result */
        for (unsigned int i = 0; i < N; ++i)
        {
            dvolt[i]->addValue(circuit()->getNodeVolt(getInterestNode(i)));
        }
        for (unsigned int i = 0; i < M; ++i)
        {
            dcurrent[i]->addValue(circuit()->getBranchCurrent(getInterestBranch(i)));
        }

        return 0;
    }

    int AnalyzerOP::prepareMNA()
    {
        for (auto &mif : circuit()->netlist()->models())
        {
            /*auto model0 = mif.model;
            std::string str(model0->name());
            testMultiLog(m_log,"model0->name():"+str);
            testMultiLog(m_log,"model0->getNumTerml():"+std::to_string(model0->getNumTerml()));
            testMultiLog(m_log,"model0->getNumBranches():"+std::to_string(model0->getNumBranches()));
            testMultiLog(m_log,"model0->getNumInnerNode():"+std::to_string(model0->getNumInnerNode()));
            for(auto &a:model0->get_m_nodes())
            {
                testMultiLog(m_log,"prepareMNA-get_m_nodes,"+std::to_string(a));
            }
           for(auto &a:model0->get_m_branches())
            {
                testMultiLog(m_log,"prepareMNA-get_m_branches,"+std::to_string(a));
            }
            auto entry0 = mif.entry;*/
            UPDATE_RC(mif.model->prepareOP());
        }
        return 0;
    }

    int AnalyzerOP::iterateMNA()
    {
        int rc = 0;
        for (auto &mif : circuit()->netlist()->models())
        {
            int ret = mif.model->iterateOP();
            if (ret == CERR_NON_CONVERGENCE)
                rc = ret;
            else if (CSIM_FAILED(ret))
                return ret;
        }
        return rc;
    }

}
