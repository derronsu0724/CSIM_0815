#include <stdio.h>
#include "constants.h"
#include "csim/utils/errors.h"
#include "csim/utils/constants.h"
#include "csim/internal/Analyzers.h"
#include "csim/internal/Circuit.h"
#include "csim/internal/Dataset.h"
#include "csim/internal/Netlist.h"
#include "csim/internal/ModelLoader.h"
#include <algorithm>
#include <complex>
#include <cstring>
#include <iostream>  
#include <fstream>  
#include <string>  

#include "log.h"
#include "parse.h"
/*
// 函数，找出每个节点连接的边
std::map<std::string, std::vector<std::pair<spice::Edge, std::string>>> findEdgesForNodes(const std::vector<spice::Edge>& edges) {  
    std::map<std::string, std::vector<std::pair<spice::Edge, std::string>>> nodeEdges;
    for (const auto& edge : edges) {
        if (!edge.type.empty()) {
            // 添加到 from 连接
            nodeEdges[edge.from].emplace_back(edge, "0");
            // 添加到 to 连接
            nodeEdges[edge.to].emplace_back(edge, "1");
        }
    }
    // 输出结果  
    
    for (const auto& node : nodeEdges) {  
        std::cout << "Node: " << node.first << "\n";
        for (const auto& edgeInfo : node.second) {
                const spice::Edge& edge = edgeInfo.first;
                const std::string& position = edgeInfo.second;
                std::cout << "  Edge: " << edge.from << " -> " << edge.to
                        << " (Type: " << edge.type << ", component: " << edge.component
                        << ") " << position << "\n";
        }
    }
    
    return nodeEdges;  
} 

// 函数，生成每个节点下所有边的两两配对  
void printEdgePairs(const std::map<std::string, std::vector<std::pair<spice::Edge, std::string>>>& nodeEdges) {  
    for (const auto& node : nodeEdges) {  
        std::cout << "Node: " << node.first << "\n";
        const auto& edgesInfo = node.second;
        size_t count = edgesInfo.size();        
        for (size_t i = 0; i < count; ++i) {  
            for (size_t j = i + 1; j < count; ++j) {  
                // 两条边的配对  
                const spice::Edge& edge1 = edgesInfo[i].first;  
                const spice::Edge& edge2 = edgesInfo[j].first;
                const std::string& position1 = edgesInfo[i].second;
                const std::string& position2 = edgesInfo[j].second;
                std::cout << edge1.component << ", "<< position1 << ", " << edge2.component << ", " << position2<< "\n";
                
                std::cout << "  Pair: ("   
                          << edge1.from << " -> " << edge1.to   
                          << ") and ("   
                          << edge2.from << " -> " << edge2.to   
                          << ") (Components: "   
                          << edge1.component << ", " << edge2.component   
                          << ", Types: " << edge1.type << ", " << edge2.type   
                          << ", Values: " << edge1.value << ", " << edge2.value
                          << ", position1: " << position1 << "," <<"position2: " << position2
                          << ")\n";
            }
        }
    }
}
std::vector<std::pair<std::string, int>> findEdgesBetweenProbes(std::vector<spice::Edge> edges,std::vector<std::string> probe_names)
{
    std::vector<std::pair<std::string, int>> componentFromPair(probe_names.size());
    int ii=0;
    for (const auto& probe_1  : probe_names) {
        // std::cout << "Probe name: " << probe_1 << ",ii: " << ii<< std::endl;  
        for (const auto& edge : edges) {
            if ( edge.from == probe_1  )
            {
                componentFromPair[ii] = {edge.component, 0};
                break;
            } else if ( edge.to == probe_1  )
            {  
                componentFromPair[ii] = {edge.component, 1};
                break;
            }
         
        }
        ii=ii+1;
    }
    
    for (const auto& ii  : componentFromPair) {
                std::cout << "Component: " << ii.first << std::endl;  
                std::cout << "node: " << ii.second << std::endl;
    }
    return componentFromPair;
}
*/
static void OPCircuit_helper(std::set<std::string> nodes,std::vector<spice::Edge> edges,std::vector<std::string> probe_names)
{
        int ret = 0;
        csim::ModelEntry *e_R = csim::ModelLoader::load(resistorLibrary);
        csim::ModelEntry *e_VDC = csim::ModelLoader::load(VDCLibrary);
        csim::ModelEntry *e_L = csim::ModelLoader::load(InductorLibrary);
        csim::ModelEntry *e_CAP = csim::ModelLoader::load(CapacitorLibrary);
        csim::Circuit *circuit = new csim::Circuit();
        // std::cout << "\n edge:" << std::endl;
        for (auto &edge : edges) {
            if (edge.type.empty()) {
                // std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (value: " << edge.value << ")" << std::endl;  
            } else if (edge.type == "voltage") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value1: " << edge.value["DC"]<< ", value2: " << edge.value["AC"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_VDC);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "V", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(edge.value["DC"]));
            } else if (edge.type == "resistor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_R);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(edge.value["value"]));
            } else if (edge.type == "capacitor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_CAP);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "C", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(edge.value["value"]));
            } else if (edge.type == "inductor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_L);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "L", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(edge.value["value"]));
            }
        }        
        // 调用函数  
        std::map<std::string, std::vector<std::pair<spice::Edge, std::string>>> nodeEdges = spice::findEdgesForNodes(edges);
        // printEdgePairs(nodeEdges);
        ret = circuit->netlist()->prepare();
        for (const auto& node : nodeEdges) {  
            //std::cout << "Node: " << node.first << "\n";
            const auto& edgesInfo = node.second;
            size_t count = edgesInfo.size();        
            for (size_t i = 0; i < count; ++i) {  
                for (size_t j = i + 1; j < count; ++j) {  
                    // 两条边的配对  
                    const spice::Edge& edge1 = edgesInfo[i].first;  
                    const spice::Edge& edge2 = edgesInfo[j].first;
                    const std::string& position1 = edgesInfo[i].second;
                    const std::string& position2 = edgesInfo[j].second;
                    //std::cout << edge1.component << ", "<< position1 << ", " << edge2.component << ", " << position2<< "\n";                            
                    ret = circuit->netlist()->wire(edge1.component.c_str(), std::stoi(position1), edge2.component.c_str(), std::stoi(position2));            
                }
            }
        }
        ret = circuit->netlist()->generateNodes();
        /* DC analysis */
        csim::AnalyzerBase *analyzer = csim::Analyzers::createInstance("OP", circuit);
        csim::Dataset dset;
        ret = analyzer->analyze(&dset);
        /* Get nodes */
        auto probe_names_data=spice::findEdgesBetweenProbes(edges,probe_names);
        unsigned int n_gnd;
        unsigned int n3[probe_names_data.size()-1];
        ret = circuit->netlist()->getTermlNode(probe_names_data[probe_names_data.size()-1].first.c_str(), probe_names_data[probe_names_data.size()-1].second, &n_gnd);
        for (size_t ii=0;ii<probe_names_data.size()-1;ii++)
        {
            unsigned int n1;
            //std::cout << "Component: " << probe_names_data[ii].first << std::endl;
            //std::cout << "node: " << probe_names_data[ii].second << std::endl;
            ret = circuit->netlist()->getTermlNode(probe_names_data[ii].first.c_str(), probe_names_data[ii].second, &n3[ii]);
            csimModel::MComplex volt = circuit->getNodeVolt(n3[ii]) - circuit->getNodeVolt(n_gnd);
            std::cout  <<"volt"<<ii<<":" << std::abs(volt)  <<"\n";
        }
        /* Check status of Circuit object */
        // EXPECT_LT(std::abs(csimModel::MComplex(4.0, 0) - volt), epsilon_linear);
        delete circuit;
        delete e_R;
        delete e_VDC;
}

static int ACLinearCircuit(std::set<std::string> nodes,std::vector<spice::Edge> edges,std::vector<std::string> probe_names,const char *fspace, double fstart, double fstop)
{
    int ret = 0;
    csim::ModelEntry *e_R = csim::ModelLoader::load(resistorLibrary);
    csim::ModelEntry *e_L = csim::ModelLoader::load(InductorLibrary);
    csim::ModelEntry *e_CAP = csim::ModelLoader::load(CapacitorLibrary);
    csim::ModelEntry *e_VAC = csim::ModelLoader::load(VACLibrary);
    csim::Circuit *circuit = new csim::Circuit();
        for (auto &edge : edges) {
            if (edge.type.empty()) {
                // std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (value: " << edge.value << ")" << std::endl;  
            } else if (edge.type == "voltage") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value1: " << edge.value["DC"]<< ", value2: " << edge.value["AC"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_VAC);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "Vp", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(edge.value["AC"]));
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "freq", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(50.0));
            } else if (edge.type == "resistor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_R);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(edge.value["value"]));
            } else if (edge.type == "capacitor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_CAP);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "C", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(edge.value["value"]));
            } else if (edge.type == "inductor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_L);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "L", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(edge.value["value"]));
            }
        } 
        // 调用函数  
        std::map<std::string, std::vector<std::pair<spice::Edge, std::string>>> nodeEdges = spice::findEdgesForNodes(edges);
        // printEdgePairs(nodeEdges);
        ret = circuit->netlist()->prepare();
        for (const auto& node : nodeEdges) {  
           // std::cout << "Node: " << node.first << "\n";
            const auto& edgesInfo = node.second;
            size_t count = edgesInfo.size();        
            for (size_t i = 0; i < count; ++i) {  
                for (size_t j = i + 1; j < count; ++j) {  
                    // 两条边的配对  
                    const spice::Edge& edge1 = edgesInfo[i].first;  
                    const spice::Edge& edge2 = edgesInfo[j].first;
                    const std::string& position1 = edgesInfo[i].second;
                    const std::string& position2 = edgesInfo[j].second;
                    //std::cout << edge1.component << ", "<< position1 << ", " << edge2.component << ", " << position2<< "\n";                            
                    ret = circuit->netlist()->wire(edge1.component.c_str(), std::stoi(position1), edge2.component.c_str(), std::stoi(position2));            
                }
            }
        }
        ret = circuit->netlist()->generateNodes();

        /* AC analysis */
        csim::AnalyzerBase *analyzer = csim::Analyzers::createInstance("AC", circuit);
        analyzer->property().setProperty("fstart", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(fstart));
        analyzer->property().setProperty("fstop", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(fstop));
        analyzer->property().setProperty("fpoints", csimModel::Variant(csimModel::Variant::VariantUint32).setUint32(50));
        analyzer->property().setProperty("fspace", csimModel::Variant(csimModel::Variant::VariantString).setString(fspace));

        /* Get nodes */
        auto probe_names_data=spice::findEdgesBetweenProbes(edges,probe_names);
        unsigned int n_gnd, n1;
        ret = circuit->netlist()->getTermlNode(probe_names_data[0].first.c_str(), probe_names_data[0].second, &n1);
        ret = circuit->netlist()->getTermlNode(probe_names_data[probe_names_data.size()-1].first.c_str(), probe_names_data[probe_names_data.size()-1].second, &n_gnd);
        analyzer->addInterestNode(n_gnd);
        analyzer->addInterestNode(n1);

        csim::Dataset dset;
        ret = analyzer->analyze(&dset);
        /* Check solution vector of AC analyzer */
        csim::Variable &F = dset.getIndependentVar("frequency");
        csim::Variable &Vgnd = dset.getDependentVar("voltage", analyzer->makeVarName("V", n_gnd));
        csim::Variable &Vn1 = dset.getDependentVar("voltage", analyzer->makeVarName("V", n1));
 
        std::ofstream outFile("output.csv"); // 创建一个 ofstream 对象并打开文件 
        for(size_t i=0; i<F.getNumValues(); ++i) {

            csimModel::MComplex _volt = Vn1.at(i) - Vgnd.at(i);
            auto volt = std::complex(_volt.real(), _volt.imag());
            outFile  <<i<<"," <<F.at(i).real()<<","<<std::abs(volt) <<"\n";
            double omega = F.at(i).real() * 2 * M_PI;
            //std::complex<double> z(R, omega*L - 1.0/(omega * C));
            //auto Ic = Vp / z;
            //double mag = std::abs(R*Ic);
            //XPECT_NEAR(mag, std::abs(volt), epsilon_linear);
            //EXPECT_NEAR(std::arg(volt), std::arg(Ic), epsilon_linear);
        }
        outFile.close(); // 关闭文件
        delete analyzer;
        delete circuit;
        delete e_R;
        delete e_L;
        delete e_CAP;
        delete e_VAC;
    return ret;
}



int main(int argc, char *argv[]) {
    std::vector<std::string> temp1(argc);
    for (int i = 0; i < argc; i++) {
        temp1.at(i) = argv[i];
        //std::cout<<temp1.at(i)<<"\n";
    }
    spice::parse_ parser(temp1.at(1));
    std::set<std::string> nodes=parser.m_nodes;
    std::vector<spice::Edge> edges=parser.m_edges;
    std::vector<std::string> probe_names=parser.m_probe_names;
    std::unordered_multimap<std::string, std::vector<double>>  analysis_value=parser.m_analysis;
    if (argc > 2) {
        std::string option = argv[2]; // 获取第一个参数（参数索引从1开始）  
        if (option == "--verbose") 
        {    // 使用范围for循环打印  
            std::cout << "Nodes in set:" << std::endl;  
            for (const auto& node : nodes) {  
                std::cout << node << std::endl;  
            }
            std::cout << "Probe :" << std::endl;  
            for (const auto& node : probe_names) {  
                std::cout << node << std::endl;  
            }
                for (auto &edge : edges) {
                    if (edge.type.empty()) {
                        // std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (value: " << edge.value << ")" << std::endl;  
                    } else if (edge.type == "voltage") {            
                        std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value1: " << edge.value["DC"]<< ", value2: " << edge.value["AC"] << ")" << std::endl;
                    } else if (edge.type == "resistor") {            
                        std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                    } else if (edge.type == "capacitor") {            
                        std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                    } else if (edge.type == "inductor") {            
                        std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                    }
                }

            auto range = analysis_value.equal_range("ac");
            std::vector<double> freq;
            for (auto it = range.first; it != range.second; ++it) {  
            std::cout << "Key: " << it->first << ", Values: ";  
            for (const auto& v : it->second) {
                freq.push_back(v);  
                //std::cout << v << " ";  
            }  
            std::cout << std::endl;}
        }
    }
    //OPCircuit_helper(nodes,edges, probe_names);
    ACLinearCircuit(nodes,edges, probe_names,"lin", 70000, 170000);

  return 1;
}
// cmake .. -Denable_testcases=ON -Denable_coverage=ON -DCMAKE_BUILD_TYPE=Debug -G "MSYS Makefiles"
// make