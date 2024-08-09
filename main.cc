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
    std::string value;
};

struct Subcircuit {
    std::string name; //  子电路名  
    std::vector<std::string> nodes;
    std::vector<std::string> external_nodes;
    std::vector<Edge> components; //  子电路内部元件  
};

void fun1(std::string& line,std::vector<std::string>& node,std::string& type, std::string& value)
{
          std::istringstream iss(line);
          std::string component;
          iss >> component; // 读取元件名称
          std::string node1,node2;
          if (component[0] == 'R' ) {              
              iss >> node1 >> node2 >> value; // 读取节点和元件的值 
              type="resistor";
          } else if (component[0] == 'C') { 
              iss >> node1 >> node2 >> value; // 读取节点和元件的值  
              type="capacitor";
          } else if (component[0] == 'L') {
              iss >> node1 >> node2 >> value; // 读取节点和元件的值  
              type="inductor";
          } else if (component[0] == 'V') { // 处理电压源 
              iss >> node1 >> node2 >> type >> value; // 读取节点、类型和电压值  
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
                        << " (Type: " << edge.type << ", Value: " << edge.value<< ", component: " << edge.component
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

static void DCCircuit_helper(std::set<std::string> nodes,std::vector<Edge> edges)
{
        int ret = 0;
        csim::ModelEntry *e_R = csim::ModelLoader::load(resistorLibrary);
        csim::ModelEntry *e_VDC = csim::ModelLoader::load(VDCLibrary);
        csim::ModelEntry *e_L = csim::ModelLoader::load(InductorLibrary);
        csim::ModelEntry *e_CAP = csim::ModelLoader::load(CapacitorLibrary);
        csim::Circuit *circuit = new csim::Circuit();
        // 打印节点
        std::cout << "node:" ;
        for (const auto &node : nodes) {  
            std::cout << node << ",";
        }
        std::cout << std::endl;
        // 打印边
        std::cout << "\n edge:" << std::endl;
        for (const auto &edge : edges) {
            if (edge.type.empty()) {
                std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (value: " << edge.value << ")" << std::endl;  
            } else if (edge.type == "DC") {            
                std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_VDC);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "V", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::stod(edge.value)));
            } else if (edge.type == "resistor") {            
                std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_R);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::stod(edge.value)));
            } else if (edge.type == "capacitor") {            
                std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_CAP);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "C", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::stod(edge.value)));
            } else if (edge.type == "inductor") {            
                std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_L);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "L", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::stod(edge.value)));
            }
        }
        // 调用函数  
        std::map<std::string, std::vector<std::pair<Edge, std::string>>> nodeEdges = findEdgesForNodes(edges);
        printEdgePairs(nodeEdges);
}

int main(int argc, char *argv[]) {
    std::vector<std::string> temp1(argc);
    for (int i = 0; i < argc; i++) {
        temp1.at(i) = argv[i];
    }
    if(temp1.at(1) == "1")
    {
        std::cout   <<__LINE__  <<"\n";
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
        int ret = 0;
        csim::ModelEntry *e_R = csim::ModelLoader::load(resistorLibrary);
        csim::ModelEntry *e_VDC = csim::ModelLoader::load(VDCLibrary);

        csim::Circuit *circuit = new csim::Circuit();
        /* Add models */
        ret = circuit->netlist()->addComponent("R1", e_R);
        ret = circuit->netlist()->addComponent("R2", e_R);
        ret = circuit->netlist()->addComponent("V1", e_VDC);

        /* Configure */
        ret = circuit->netlist()->configComponent("R1", "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(20.0));
        ret = circuit->netlist()->configComponent("R2", "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(10.0));
        ret = circuit->netlist()->configComponent("V1", "V", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(12.0));

        /* Set netlist */
        ret = circuit->netlist()->prepare();
        ret = circuit->netlist()->wire("R1", 0, "V1", 0);
        ret = circuit->netlist()->wire("R1", 1, "R2", 0);
        ret = circuit->netlist()->wire("R2", 1, "V1", 1);
        ret = circuit->netlist()->generateNodes();

        /* DC analysis */
        csim::AnalyzerBase *analyzer = csim::Analyzers::createInstance("OP", circuit);
        csim::Dataset dset;
        ret = analyzer->analyze(&dset);
        /* Get nodes */
        unsigned int n_gnd, n1;
        ret = circuit->netlist()->getTermlNode("R1", 1, &n1);
        ret = circuit->netlist()->getTermlNode("R2", 1, &n_gnd);
        /* Check status of Circuit object */
        csimModel::MComplex volt = circuit->getNodeVolt(n1) - circuit->getNodeVolt(n_gnd);
        std::cout  <<"volt:" << std::abs(volt)  <<"\n";
        // EXPECT_LT(std::abs(csimModel::MComplex(4.0, 0) - volt), epsilon_linear);
        delete circuit;
        delete e_R;
        delete e_VDC;

    } else if  (temp1.at(1) == "2") {
        Graph graph(6);      
        graph.addEdge(0, 1);
        graph.addEdge(0, 2);
        graph.addEdge(1, 3);
        graph.addEdge(1, 4);
        graph.addEdge(2, 5);      
        std::cout << "DFS traversal: ";
        graph.DFS(0);
    }
    else if  (temp1.at(1) == "3") {
        std::ofstream netlist("circuit.sp");
        // 检查文件是否成功打开  
        if (!netlist) {  
            std::cerr << "无法打开文件！" << std::endl;
            return 1;
        }
        // 写入网表内容  
        netlist << "* Sample HSPICE Netlist\n";  
        netlist << "V1 1 0 DC 5V\n";           // 电压源  
        netlist << "R1 1 2 1k\n";              // 电阻  
        netlist << "C1 2 0 10n\n";             // 电容  
        netlist << "L2 2 0 11\n"; 
        // netlist << ".model NMOS NMOS (VTO=0.7)\n"; // NMOS模型定义  
         // netlist << "M1 3 1 0 0 NMOS W=10u L=1u\n"; // 晶体管  
        netlist << ".tran 0.1ms 10ms\n";       // 转换分析指令  
        netlist << ".end\n";                   // 网表结束
        // 关闭文件  
        netlist.close();
        std::cout << "finish" << std::endl;
    }
    else if  (temp1.at(1) == "4") {
        std::ifstream netlist("circuit.sp");  
        std::string line; 
        std::vector<Subcircuit> subcircuits; // 存储所有子电路
        Subcircuit current_subcircuit;  
        bool in_subcircuit = false;
        std::set<std::string> nodes;  // 存储节点  
        std::vector<Edge> edges;       // 存储边
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
            std::string value; // 存储元件的值
            // 解析节点和元件值
            if (component == ".SUBCKT") { // 检测子电路开始
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
        /*
        // 打印子电路信息  
        std::cout << "\nSubcircuits:" << std::endl;  
        for (const auto &subcircuit : subcircuits) {  
            std::cout << "Subcircuit: " << subcircuit.name << std::endl;  
            std::cout << "External Nodes: ";        
            for (const auto &node : subcircuit.external_nodes) {  
                std::cout << node << " ";  
            }  
            std::cout << std::endl;  

            std::cout << "Components: " << std::endl;  
            for (const auto &component : subcircuit.components) {  
                std::cout << component.component << ": " << component.from << " -> " << component.to  
                        << " (type: " << component.type << ", value: " << component.value << ")" << std::endl;  
            }  
        }
        */
        DCCircuit_helper(nodes,edges);
    }
  return 1;
}
// cmake .. -Denable_testcases=ON -Denable_coverage=ON -DCMAKE_BUILD_TYPE=Debug -G "MSYS Makefiles"
// make