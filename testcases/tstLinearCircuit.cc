#include "gtest/gtest.h"
#include "constants.h"
#include "csim/utils/errors.h"
#include "csim/internal/Circuit.h"
#include "csim/internal/Netlist.h"
#include "csim/internal/ModelLoader.h"
#include <algorithm>
#include <cstring>

namespace csim
{
    /*
     * Circuit diagram
     *         20ohm           10ohm
     *     0 .______. 1    0 .______. 1
     *  +---+|__R1__|+------+|__R2__|+---+
     *  |                V1              |
     *  |          0    ,-.    1         |
     *  +--------------(---)-------------+
     *                  `-'
     *                  12V
     */
    TEST(circuit_R_VS_divider, tstLinearCircuit)
    {
        int ret = 0;
        ModelEntry *e_R = ModelLoader::load(resistorLibrary);
        ASSERT_NE(nullptr, e_R);
        ModelEntry *e_VDC = ModelLoader::load(VDCLibrary);
        ASSERT_NE(nullptr, e_VDC);

        Circuit *circuit = new Circuit();

        /* Add models */
        ret = circuit->netlist()->addComponent("R1", e_R);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->addComponent("R2", e_R);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->addComponent("V1", e_VDC);
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        /* Configure */
        ret = circuit->netlist()->configComponent("R1", "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(20.0));
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->configComponent("R2", "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(10.0));
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->configComponent("V1", "V", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(12.0));
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        /* Set netlist */
        ret = circuit->netlist()->prepare();
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        ret = circuit->netlist()->wire("R1", 0, "V1", 0);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R1", 1, "R2", 0);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R2", 1, "V1", 1);
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        ret = circuit->netlist()->generateNodes();
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        /* DC analysis */
        ret = circuit->analyseDC();
        EXPECT_EQ(CERR_SUCCEEDED, ret);

        /* Check voltages */
        unsigned int n_gnd, n1;
        ret = circuit->netlist()->getTermlNode("R1", 1, &n1);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("R2", 1, &n_gnd);
        EXPECT_EQ(CERR_SUCCEEDED, ret);

        Complex volt = circuit->getNodeVolt(n1) - circuit->getNodeVolt(n_gnd);
        EXPECT_LT(std::abs(Complex(4.0, 0) - volt), epsilon);

        delete circuit;
        delete e_R;
        delete e_VDC;
    }

    /*
     *                                +-----------+-----+
     *   Circuit diagram:            1|          0|     |0
     *                               .-.         .-.   .-.
     *                         30ohm |R|   40ohm |R|   |R| 60ohm
     *                               |3|         |4|   |6|
     *                               '-'         '-'   '-'
     *                               0|          1|     |1
     *          10ohm         20ohm   |    50ohm  |     |
     *      0 .______. 1  0 .______. 1| 0.______.1|     |
     *   +----|__R1__|------|__R2__|--+--|__R5__|-+-----+
     *   |                V1                            |
     *   |          0    ,-.    1                       |
     *   +--------------(---)---------------------------+
     *                   `-'
     *                   12V
     */
    TEST(circuit_R_VS_network, tstLinearCircuit)
    {
        int ret = 0;
        ModelEntry *e_R = ModelLoader::load(resistorLibrary);
        ASSERT_NE(nullptr, e_R);
        ModelEntry *e_VDC = ModelLoader::load(VDCLibrary);
        ASSERT_NE(nullptr, e_VDC);

        Circuit *circuit = new Circuit();

        /* Add 6 resistors */
        char ref[32];
        for (int i = 1; i <= 6; ++i)
        {
            snprintf(ref, sizeof(ref), "R%d", i);
            ret = circuit->netlist()->addComponent(ref, e_R);
            ASSERT_EQ(CERR_SUCCEEDED, ret);
        }
        ret = circuit->netlist()->addComponent("V1", e_VDC);
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        /* Configure */
        for (int i = 1; i <= 6; ++i)
        {
            snprintf(ref, sizeof(ref), "R%d", i);
            ret = circuit->netlist()->configComponent(ref, "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(10.0 * i));
            ASSERT_EQ(CERR_SUCCEEDED, ret);
        }
        ret = circuit->netlist()->configComponent("V1", "V", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(12.0));
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        /* Set netlist */
        ret = circuit->netlist()->prepare();
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        ret = circuit->netlist()->wire("R1", 0, "V1", 0);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R1", 1, "R2", 0);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R2", 1, "R5", 0);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R2", 1, "R3", 0);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R3", 1, "R4", 0);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R3", 1, "R6", 0);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R5", 1, "R4", 1);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R5", 1, "R6", 1);
        ASSERT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->wire("R5", 1, "V1", 1);
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        ret = circuit->netlist()->generateNodes();
        ASSERT_EQ(CERR_SUCCEEDED, ret);

        /* Test nodes */
        EXPECT_EQ(5, circuit->netlist()->getNumNodes());

        unsigned int nds[5];
        ret = circuit->netlist()->getTermlNode("R1", 0, &nds[0]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("V1", 0, &nds[1]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        EXPECT_EQ(nds[0], nds[1]);

        ret = circuit->netlist()->getTermlNode("R1", 1, &nds[0]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("R2", 0, &nds[1]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        EXPECT_EQ(nds[0], nds[1]);

        ret = circuit->netlist()->getTermlNode("R2", 1, &nds[0]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("R3", 0, &nds[1]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("R5", 0, &nds[2]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        EXPECT_EQ(nds[0], nds[1]);
        EXPECT_EQ(nds[1], nds[2]);

        ret = circuit->netlist()->getTermlNode("R3", 1, &nds[0]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("R4", 0, &nds[1]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("R6", 0, &nds[2]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        EXPECT_EQ(nds[0], nds[1]);
        EXPECT_EQ(nds[1], nds[2]);

        ret = circuit->netlist()->getTermlNode("R5", 1, &nds[0]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("R4", 1, &nds[1]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("R6", 1, &nds[2]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        ret = circuit->netlist()->getTermlNode("V1", 1, &nds[3]);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        EXPECT_EQ(nds[0], nds[1]);
        EXPECT_EQ(nds[1], nds[2]);
        EXPECT_EQ(nds[2], nds[3]);

        /* DC analysis */
        ret = circuit->analyseDC();
        EXPECT_EQ(CERR_SUCCEEDED, ret);

        /* Check voltages */
        unsigned int n_gnd, n1;
        ret = circuit->netlist()->getTermlNode("V1", 1, &n_gnd);
        EXPECT_EQ(CERR_SUCCEEDED, ret);

        Complex volt;

        ret = circuit->netlist()->getTermlNode("R1", 1, &n1);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        volt = circuit->getNodeVolt(n1) - circuit->getNodeVolt(n_gnd);
        EXPECT_LT(std::abs(Complex(9.85567, 0) - volt), epsilon);

        ret = circuit->netlist()->getTermlNode("R2", 1, &n1);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        volt = circuit->getNodeVolt(n1) - circuit->getNodeVolt(n_gnd);
        EXPECT_LT(std::abs(Complex(5.56701, 0) - volt), epsilon);

        ret = circuit->netlist()->getTermlNode("R3", 1, &n1);
        EXPECT_EQ(CERR_SUCCEEDED, ret);
        volt = circuit->getNodeVolt(n1) - circuit->getNodeVolt(n_gnd);
        EXPECT_LT(std::abs(Complex(2.47423, 0) - volt), epsilon);

        delete circuit;
        delete e_R;
        delete e_VDC;
    }

} // namespace
