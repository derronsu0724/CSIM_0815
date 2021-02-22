#ifndef RESISTOR_H_
#define RESISTOR_H_

#include "csim/model/ModelBase.h"

namespace csimModel
{

    class Resistor : public ModelBase
    {
    public:
        Resistor(MODELBASE_CONSTRUCTOR_DEF);
        virtual ~Resistor();

    public:
        virtual int configure();
        virtual int prepareDC();
        virtual int prepareAC();
        virtual int prepareTR();
        virtual int iterateDC();
        virtual int iterateAC();
        virtual int iterateTR();

    private:
        double m_G;
    };

}

#endif // RESISTOR_H_