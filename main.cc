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
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "V", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["DC"])));
            } else if (edge.type == "resistor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_R);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["value"])));
            } else if (edge.type == "capacitor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_CAP);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "C", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["value"])));
            } else if (edge.type == "inductor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_L);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "L", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["value"])));
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
        unsigned int n_gnd,n1[probe_names_data.size()-1];
        ret = circuit->netlist()->getTermlNode(probe_names_data[probe_names_data.size()-1].first.c_str(), probe_names_data[probe_names_data.size()-1].second, &n_gnd);
        for (size_t ii=0;ii<probe_names_data.size()-1;ii++)
        {
            ret = circuit->netlist()->getTermlNode(probe_names_data[ii].first.c_str(), probe_names_data[ii].second, &n1[ii]);
            csimModel::MComplex volt = circuit->getNodeVolt(n1[ii]) - circuit->getNodeVolt(n_gnd);
            std::cout  <<"volt,"<<probe_names[ii]<<":" << std::abs(volt)  <<"\n";
        }
        delete circuit;
        delete e_R;
        delete e_VDC;
}

static int ACLinearCircuit(std::set<std::string> nodes,std::vector<spice::Edge> edges,std::vector<std::string> probe_names,const char *fspace, double fstart, double fstop,int npoints)
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
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "Vp", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["AC"])));
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "freq", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(50.0));
            } else if (edge.type == "resistor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_R);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["value"])));
            } else if (edge.type == "capacitor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_CAP);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "C", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["value"])));
            } else if (edge.type == "inductor") {            
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_L);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "L", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["value"])));
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
        analyzer->property().setProperty("fpoints", csimModel::Variant(csimModel::Variant::VariantUint32).setUint32(npoints));
        analyzer->property().setProperty("fspace", csimModel::Variant(csimModel::Variant::VariantString).setString(fspace));

        /* Get nodes */
        auto probe_names_data=spice::findEdgesBetweenProbes(edges,probe_names);
        unsigned int n_gnd, n1[probe_names_data.size()-1];
        
        ret = circuit->netlist()->getTermlNode(probe_names_data[probe_names_data.size()-1].first.c_str(), probe_names_data[probe_names_data.size()-1].second, &n_gnd);
        analyzer->addInterestNode(n_gnd);
        for (size_t ii=0;ii<probe_names_data.size()-1;ii++)
        {
            ret = circuit->netlist()->getTermlNode(probe_names_data[ii].first.c_str(), probe_names_data[ii].second, &n1[ii]);
            analyzer->addInterestNode(n1[ii]);
        }
        csim::Dataset dset;
        ret = analyzer->analyze(&dset);
        /* Check solution vector of AC analyzer */
        csim::Variable &F = dset.getIndependentVar("frequency");
        csim::Variable &Vgnd = dset.getDependentVar("voltage", analyzer->makeVarName("V", n_gnd));
        std::cout<<"freq numbers: "<<F.getNumValues()<<"\n";
        std::cout<<"probe numbers: "<<probe_names_data.size()-1<<"\n";
        std::vector<std::vector<std::complex<double>>> result(F.getNumValues(), std::vector<std::complex<double>>(probe_names_data.size()-1,std::complex<double>(0,0)));
        std::ofstream outFile("output.csv"); // 创建一个 ofstream 对象并打开文件 
        for(size_t i=0; i<F.getNumValues(); ++i) {
            outFile  <<i<<","<<F.at(i).real()<<",";
            for (size_t jj=0;jj<probe_names_data.size()-1;jj++)
            {
                csim::Variable &Vn1 = dset.getDependentVar("voltage", analyzer->makeVarName("V", n1[jj]));
                csimModel::MComplex _volt = Vn1.at(i) - Vgnd.at(i);
                result[i][jj] = std::complex(_volt.real(), _volt.imag());
                auto volt = result[i][jj];
                outFile<<std::abs(volt) <<",";

            }
            outFile  <<"\n";
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

static void tstTransient(std::set<std::string> nodes,std::vector<spice::Edge> edges,std::vector<std::string> probe_names, double tstep, double tstop)
{
        //std::cout <<"tstep "<< tstep << "\n";
        //std::cout <<"tstop "<< tstop << "\n";
        int ret = 0;
        csim::ModelEntry *e_D = csim::ModelLoader::load(PNLibrary);
        csim::ModelEntry *e_R = csim::ModelLoader::load(resistorLibrary);
        csim::ModelEntry *e_L = csim::ModelLoader::load(InductorLibrary);
        csim::ModelEntry *e_CAP = csim::ModelLoader::load(CapacitorLibrary);
        csim::ModelEntry *e_VAC = csim::ModelLoader::load(VACLibrary);
        csim::Circuit *circuit = new csim::Circuit();

        for (auto &edge : edges) {
            if (edge.type.empty()) {
                // std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (value: " << edge.value << ")" << std::endl;  
            } else if (edge.type == "voltage") {
                double vp =std::get<std::vector<double>>(edge.value["sin"])[1];
                double freq=std::get<std::vector<double>>(edge.value["sin"])[2];
                double phase=std::get<std::vector<double>>(edge.value["sin"])[5];
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value1: " << edge.value["DC"]<< ", value2: " << edge.value["AC"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_VAC);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "Vp", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(vp));
                //ret = circuit->netlist()->configComponent(edge.component.c_str(), "freq", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(50.0));
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "freq", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(freq));
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "phase", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(phase));
                /*
                    // 输出解析结果
                    auto value0=edge.value;
                    auto pwl_data=std::get<std::vector<spice::TimeVoltage>>(edge.value["PWL"]);
                    for (const auto& pair : value0) {
                        std::cout << "Key: " << pair.first << "\n";
                        if (std::holds_alternative<double>(pair.second)) {
                            std::cout << "Value: " << std::get<double>(pair.second) << "\n";
                        } else if (std::holds_alternative<std::vector<spice::TimeVoltage>>(pair.second)) {
                            for (const auto& tv : std::get<std::vector<spice::TimeVoltage>>(pair.second)) {
                                std::cout << "Time: " << tv.time << ", Voltage: " << tv.voltage << "\n";
                            }
                        }
                    }
                */
            } else if (edge.type == "resistor") {
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_R);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["value"])));
            } else if (edge.type == "capacitor") {
                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_CAP);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "C", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["value"])));
            } else if (edge.type == "inductor") {

                //std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value["value"] << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_L);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "L", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::get<double>(edge.value["value"])));
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
        /* Transient analysis */
        csim::AnalyzerBase *analyzer = csim::Analyzers::createInstance("transient", circuit);
        analyzer->property().setProperty("tstop", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(tstop));
        analyzer->property().setProperty("tstep", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(tstep));

        /* Get nodes */
        auto probe_names_data=spice::findEdgesBetweenProbes(edges,probe_names);
        unsigned int n_gnd, n1[probe_names_data.size()-1];        
        ret = circuit->netlist()->getTermlNode(probe_names_data[probe_names_data.size()-1].first.c_str(), probe_names_data[probe_names_data.size()-1].second, &n_gnd);
        analyzer->addInterestNode(n_gnd);
        for (size_t ii=0;ii<probe_names_data.size()-1;ii++)
        {
            ret = circuit->netlist()->getTermlNode(probe_names_data[ii].first.c_str(), probe_names_data[ii].second, &n1[ii]);
            analyzer->addInterestNode(n1[ii]);
        }

        csim::Dataset dset;
        ret = analyzer->analyze(&dset);
        /* Check solution vector of DC analyzer */
        const csim::Variable &T = dset.getIndependentVar("time");
        csim::Variable &Vgnd = dset.getDependentVar("voltage", analyzer->makeVarName("V", n_gnd));
        //const csim::Variable &Vn1 = dset.getDependentVar("voltage", analyzer->makeVarName("V", n1));
        std::ofstream fof("outsD.csv");
        for (unsigned int i = 0; i < T.getNumValues(); ++i)
        {

            fof  <<i<<","<<T.at(i).real()<<",";
            for (size_t jj=0;jj<probe_names_data.size()-1;jj++)
            {
                csim::Variable &Vn1 = dset.getDependentVar("voltage", analyzer->makeVarName("V", n1[jj]));
                csimModel::MComplex _volt = Vn1.at(i) - Vgnd.at(i);
                fof << _volt.real() << ",";
            }
            fof  <<"\n";
        }

        delete analyzer;
        delete circuit;
        delete e_D;
        delete e_R;
        delete e_VAC;
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
                        std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value1: " << std::get<double>(edge.value["DC"])<< ", value2: " << std::get<double>(edge.value["AC"]) << ")" << std::endl;
                    } else if (edge.type == "resistor") {            
                        std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << std::get<double>(edge.value["value"]) << ")" << std::endl;
                    } else if (edge.type == "capacitor") {            
                        std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << std::get<double>(edge.value["value"])<< ")" << std::endl;
                    } else if (edge.type == "inductor") {            
                        std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << std::get<double>(edge.value["value"] )<< ")" << std::endl;
                    }
                }
            std::vector<std::string> ana{"ac","op","tran"};
            for(auto ii:ana){
                auto range = analysis_value.equal_range(ii);                
                for (auto it = range.first; it != range.second; ++it) {  
                std::cout << "Key: " << it->first << ", Values: ";  
                for (const auto& v : it->second) {
                    std::cout << v << " ";  
                }  
                std::cout << std::endl;
                }
            }
        }
    }
    for (auto it = analysis_value.begin(); it != analysis_value.end(); ++it) {
        std::cout << "Key: " << it->first << ", Value: ";
        if(it->first=="ac"){
            std::vector<double> freq;
            for (double value : it->second) {
                freq.push_back(value); 
                std::cout << value << " ";
            }
            std::cout << std::endl;
            ACLinearCircuit(nodes,edges, probe_names,"lin", freq[1], freq[2],freq[0]);
        } else if(it->first=="op"){
            OPCircuit_helper(nodes,edges, probe_names);
        } else if(it->first=="tran"){
            std::vector<double> time;
            for (double value : it->second) {
                time.push_back(value); 
                std::cout << value << " ";
            }
            std::cout << "\n";
            std::cout <<"time[0],"<<time[0]<< "\n";
            std::cout <<"time[1],"<<time[1]<< "\n";
            tstTransient(nodes,edges, probe_names,time[0], time[1]);
        }
    }   


  return 1;
}
// cmake .. -Denable_testcases=ON -Denable_coverage=ON -DCMAKE_BUILD_TYPE=Debug -G "MSYS Makefiles"
// cmake .. -DCMAKE_BUILD_TYPE=Debug -G "MSYS Makefiles" -DDEBUG_MODE=ON
// make