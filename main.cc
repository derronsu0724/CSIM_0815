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
#include <boost/regex.hpp>
#include "log.h"
std::vector<std::string> split_space(const std::string& text) {
  boost::regex ws_re("\\s+");
  std::vector<std::string> vector_(
      boost::sregex_token_iterator(text.begin(), text.end(), ws_re, -1),
      boost::sregex_token_iterator());
  return vector_;
}
class Graph {
private:
    std::vector<std::vector<int>> adjacencyList;
    std::vector<bool> visited;
 
public:
    Graph(int n) {
        adjacencyList.resize(n);
        visited.resize(n, false);
    }
 
    void addEdge(int u, int v) {
        adjacencyList[u].push_back(v);
        adjacencyList[v].push_back(u);
    }
 
    void DFS(int vertex) {
        visited[vertex] = true;
        std::cout << vertex << " ";
 
        for (int i = 0; i < adjacencyList[vertex].size(); i++) {
            int neighbor = adjacencyList[vertex][i];
            if (!visited[neighbor]) {
                DFS(neighbor);
            }
        }
    }
};

struct Edge {
    std::string from;
    std::string to;
    std::string component;
    std::string type;
    std::unordered_map<std::string, double> value; // 以名称为Key的参数值,电阻电容电感值
};

struct Subcircuit {
    std::string name; //  子电路名  
    std::vector<std::string> nodes;
    std::vector<std::string> external_nodes;
    std::vector<Edge> components; //  子电路内部元件  
};

void fun1(std::string& line,std::vector<std::string>& node,std::string& type, std::unordered_map<std::string, double>& data1)
{
          std::istringstream iss(line);
          std::string component;
          iss >> component; // 读取元件名称
          std::string node1,node2;
          if (std::tolower(component[0]) == 'r' ) {
              double value;
              iss >> node1 >> node2 >> value; // 读取节点和元件的值 
              type="resistor";
              //data1.push_back(value);
              data1["value"]=value;
          } else if (std::tolower(component[0] )== 'c') {
                double value; 
              iss >> node1 >> node2 >> value; // 读取节点和元件的值
              type="capacitor";
              data1["value"]=value;
              //data1.push_back(value);
          } else if (std::tolower(component[0] )== 'l') {
              double value;
              iss >> node1 >> node2 >> value; // 读取节点和元件的值
              type="inductor";
              //data1.push_back(value);
              data1["value"]=value;
          } else if (std::tolower(component[0] )== 'v') { // 处理电压源 
              auto a1=split_space(line);
              double value1, value2;
              node1=a1[1];
              node2=a1[2];
              value1=std::stod(a1[4]);
              value2=std::stod(a1[6]);              
              type="voltage";
              //data1.push_back(value1);
              //data1.push_back(value2);
              data1["DC"]=value1;
              data1["AC"]=value2;
          }
          node.push_back(node1);
          node.push_back(node2);
}

int findNodePosition(const std::string& node1, const std::vector<std::string>& a1) {  
    for (size_t i = 0; i < a1.size(); ++i) {  
        if (a1[i] == node1) {  
            return static_cast<int>(i); // 找到时返回索引  
        }  
    }  
    return -1; // -1 表示未找到  
} 
  
// 函数，找出每个节点连接的边
std::map<std::string, std::vector<std::pair<Edge, std::string>>> findEdgesForNodes(const std::vector<Edge>& edges) {  
    std::map<std::string, std::vector<std::pair<Edge, std::string>>> nodeEdges;
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
                const Edge& edge = edgeInfo.first;
                const std::string& position = edgeInfo.second;
                std::cout << "  Edge: " << edge.from << " -> " << edge.to
                        << " (Type: " << edge.type << ", component: " << edge.component
                        << ") " << position << "\n";
        }
    }
    return nodeEdges;  
} 

// 函数，生成每个节点下所有边的两两配对  
void printEdgePairs(const std::map<std::string, std::vector<std::pair<Edge, std::string>>>& nodeEdges) {  
    for (const auto& node : nodeEdges) {  
        std::cout << "Node: " << node.first << "\n";
        const auto& edgesInfo = node.second;
        size_t count = edgesInfo.size();        
        for (size_t i = 0; i < count; ++i) {  
            for (size_t j = i + 1; j < count; ++j) {  
                // 两条边的配对  
                const Edge& edge1 = edgesInfo[i].first;  
                const Edge& edge2 = edgesInfo[j].first;
                const std::string& position1 = edgesInfo[i].second;
                const std::string& position2 = edgesInfo[j].second;
                std::cout << edge1.component << ", "<< position1 << ", " << edge2.component << ", " << position2<< "\n";
                /*
                std::cout << "  Pair: ("   
                          << edge1.from << " -> " << edge1.to   
                          << ") and ("   
                          << edge2.from << " -> " << edge2.to   
                          << ") (Components: "   
                          << edge1.component << ", " << edge2.component   
                          << ", Types: " << edge1.type << ", " << edge2.type   
                          << ", Values: " << edge1.value << ", " << edge2.value
                          << ", position1: " << position1 << "," <<"position2: " << position2
                          << ")\n";*/
            }
        }
    }
}
std::vector<std::pair<std::string, int>> findEdgesBetweenProbes(std::vector<Edge> edges,std::vector<std::string> probe_names)
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

static void OPCircuit_helper(std::set<std::string> nodes,std::vector<Edge> edges,std::vector<std::string> probe_names)
{
        probe_names.push_back("0");
        int ret = 0;
        csim::ModelEntry *e_R = csim::ModelLoader::load(resistorLibrary);
        csim::ModelEntry *e_VDC = csim::ModelLoader::load(VDCLibrary);
        csim::ModelEntry *e_L = csim::ModelLoader::load(InductorLibrary);
        csim::ModelEntry *e_CAP = csim::ModelLoader::load(CapacitorLibrary);
        csim::Circuit *circuit = new csim::Circuit();
        // 打印节点
        /*
        std::cout << "node:" ;
        for (auto &node : nodes) {  
            std::cout << node << ",";
        }
        std::cout << std::endl;*/
        // 打印边
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
        std::map<std::string, std::vector<std::pair<Edge, std::string>>> nodeEdges = findEdgesForNodes(edges);
        // printEdgePairs(nodeEdges);
        ret = circuit->netlist()->prepare();
        for (const auto& node : nodeEdges) {  
            std::cout << "Node: " << node.first << "\n";
            const auto& edgesInfo = node.second;
            size_t count = edgesInfo.size();        
            for (size_t i = 0; i < count; ++i) {  
                for (size_t j = i + 1; j < count; ++j) {  
                    // 两条边的配对  
                    const Edge& edge1 = edgesInfo[i].first;  
                    const Edge& edge2 = edgesInfo[j].first;
                    const std::string& position1 = edgesInfo[i].second;
                    const std::string& position2 = edgesInfo[j].second;
                    std::cout << edge1.component << ", "<< position1 << ", " << edge2.component << ", " << position2<< "\n";                            
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
        auto probe_names_data=findEdgesBetweenProbes(edges,probe_names);
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

static int ACLinearCircuit(std::set<std::string> nodes,std::vector<Edge> edges,std::vector<std::string> probe_names,const char *fspace, double fstart, double fstop)
{
    std::cout   <<__LINE__ <<"************ ACLinearCircuit ************\n";
    std::cout   <<fspace <<"\n";
    std::cout   <<fstart <<"\n";
    std::cout   <<fstop <<"\n";
    csim::ModelEntry *e_R = csim::ModelLoader::load(resistorLibrary);
    csim::ModelEntry *e_L = csim::ModelLoader::load(InductorLibrary);
    csim::ModelEntry *e_CAP = csim::ModelLoader::load(CapacitorLibrary);
    csim::ModelEntry *e_VAC = csim::ModelLoader::load(VACLibrary);
    csim::Circuit *circuit = new csim::Circuit();
    
    int ret = 0;
    return ret;
}



int main(int argc, char *argv[]) {
    std::vector<std::string> temp1(argc);
    for (int i = 0; i < argc; i++) {
        temp1.at(i) = argv[i];
        std::cout << temp1.at(i) << std::endl;
    }
    if(argc > 1)  {
        std::cout << temp1.at(1) << std::endl;
        std::ifstream netlist(temp1.at(1));  
        std::string line; 
        std::vector<Subcircuit> subcircuits; // 存储所有子电路
        Subcircuit current_subcircuit;  
        bool in_subcircuit = false;
        std::set<std::string> nodes;  // 存储节点  
        std::vector<Edge> edges;       // 存储边
        std::vector<std::string> probe_names;
        if (!netlist) {  
            std::cerr << "Unable to open file!" << std::endl;  
            return 1;  
        }  

        while (std::getline(netlist, line)) {
            // 跳过注释行  
            if (line.empty() || line[0] == '*') continue;
            std::istringstream iss(line);
            std::string component;
            iss >> component; // 读取元件名称
            std::vector<std::string> node_component;
            std::string type;  // 存储类型  
            //std::vector<double> value; // 存储元件的值
            std::unordered_map<std::string, double> value; // 存储元件的值
            // 解析节点和元件值
            if (component == ".subckt") { // 检测子电路开始
                in_subcircuit = true;
                iss >> current_subcircuit.name;
                std::string node;
                while (iss >> node) {
                    current_subcircuit.external_nodes.push_back(node);
                }
            } else if (in_subcircuit && component == ".ENDS") { // 子电路结束  
                subcircuits.push_back(current_subcircuit);
                current_subcircuit = Subcircuit(); // 清空当前子电路
                in_subcircuit = false;
            } else if (in_subcircuit) {// 子电路内部元件  
                fun1(line,node_component,type,value);
                current_subcircuit.components.push_back({node_component[0], node_component[1], component, type, value});  
            }
            else if (component == ".probe") {                
                // 查找 .probe 的开始位置
                size_t probe_pos = line.find(".probe");
                // 确保找到了 .probe
                if (probe_pos != std::string::npos) {
                    // 从 .probe 后面开始查找 v(
                    size_t start_pos = line.find("v(", probe_pos);
                    while (start_pos != std::string::npos) {
                        size_t end_pos = line.find(")", start_pos);
                        if (end_pos != std::string::npos) {
                            std::string variable = line.substr(start_pos + 2, end_pos - start_pos - 2);
                            probe_names.push_back(variable); // 将找到的变量保存到 vector 中
                            start_pos = line.find("v(", end_pos);
                        } else {
                            break;
                        }
                    }
                }
            }
            else if (component[0] == 'X') { 
                //std::cout   <<__LINE__ <<"\n";
                auto temp_=split_space(line);
                std::string lastElement = temp_.back();
                std::vector<std::string> a1;
                for(size_t ii=1;ii<temp_.size()-1;ii++)
                {
                    a1.push_back(temp_[ii]);
                    nodes.insert(temp_[ii]);
                }
                for (const auto &subcircuit : subcircuits) { 
                    if(subcircuit.name == lastElement)
                    {
                        for (const auto &component : subcircuit.components) {
                            std::string temp_node1,temp_node2;
                            if(findNodePosition(component.from,subcircuit.external_nodes) !=-1){
                                temp_node1=a1[findNodePosition(component.from,subcircuit.external_nodes)];
                            } else{
                                temp_node1=component.from;
                            }
                            if(findNodePosition(component.to,subcircuit.external_nodes) !=-1){
                                temp_node2=a1[findNodePosition(component.to,subcircuit.external_nodes)];
                            } else{
                                temp_node2=component.to;
                            }
                            edges.push_back({temp_node1, temp_node2, component.component+"_"+subcircuit.name, component.type, component.value});
                        }
                        break;
                    }
                    
                }
            }
            else{
            fun1(line,node_component,type,value);
            edges.push_back({node_component[0], node_component[1], component, type, value});
            nodes.insert(node_component[0]);
            nodes.insert(node_component[1]);
            }
        }  
        netlist.close();  
        OPCircuit_helper(nodes,edges, probe_names);
        //ACLinearCircuit(nodes,edges, probe_names,"lin", 500, 800);
    }
  return 1;
}
// cmake .. -Denable_testcases=ON -Denable_coverage=ON -DCMAKE_BUILD_TYPE=Debug -G "MSYS Makefiles"
// make