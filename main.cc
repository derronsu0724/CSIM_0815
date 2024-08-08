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
    std::string type; // 存储元件的类型（如 DC）  
    std::string value; // 存储元件的值（如 5V）  
};  

struct Subcircuit {  
    std::string name;      // 子电路名  
    std::vector<std::string> nodes; // 端口列表  
    std::vector<Edge> components; // 子电路内部元件  
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
        std::cout   << __LINE__  <<"\n";
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
        std::cerr << "Unable to open file！" << std::endl;  
        return 1;  
    }  

    while (std::getline(netlist, line)) {  
        // 跳过注释行  
        if (line.empty() || line[0] == '*') continue;  

        std::istringstream iss(line);  
        std::string component;
        iss >> component; // 读取元件名称
        
        std::string node1, node2;
        std::vector<std::string> node_component;
        std::string type;  // 存储类型  
        std::string value; // 存储元件的值  

        // 解析节点和元件值
        if (component == ".SUBCKT") { // 检测子电路开始  
              in_subcircuit = true;  
              iss >> current_subcircuit.name;  
              std::string node;  
              while (iss >> node) {  
                  current_subcircuit.nodes.push_back(node);  
              }  
          } else if (in_subcircuit && component == ".ENDS") { // 子电路结束  
              subcircuits.push_back(current_subcircuit);  
              current_subcircuit = Subcircuit(); // 清空当前子电路  
              in_subcircuit = false;  
          } else if (in_subcircuit) { // 子电路内部元件  
              // 在这里继续解析元件  
              fun1(line,node_component,type,value);
              current_subcircuit.components.push_back({node_component[0], node_component[1], component, type, value});  
              current_subcircuit.nodes.push_back(node_component[0]);  
              current_subcircuit.nodes.push_back(node_component[1]);  
        }
        else if (component[0] == 'X') { 
              std::cout   <<__LINE__ <<"\n";
              auto temp_=split_space(line);
              std::string lastElement = temp_.back();
              std::cout  <<"lastElement:" <<lastElement <<"\n";
              for(auto ii:temp_)
              {
                std::cout   <<ii <<"\n";
              }
        }
        else{
          fun1(line,node_component,type,value);
          edges.push_back({node_component[0], node_component[1], component, type, value});
          // 存储节点  
          nodes.insert(node_component[0]);
          nodes.insert(node_component[1]);
        }
    }  

    netlist.close();  

    // 打印子电路信息  
    std::cout << "\nSubcircuits:" << std::endl;  
    for (const auto &subcircuit : subcircuits) {  
        std::cout << "Subcircuit: " << subcircuit.name << std::endl;  
        std::cout << "Nodes: ";  
        for (const auto &node : subcircuit.nodes) {  
            std::cout << node << " ";  
        }  
        std::cout << std::endl;  

        std::cout << "Components: " << std::endl;  
        for (const auto &component : subcircuit.components) {  
            std::cout << component.component << ": " << component.from << " -> " << component.to  
                      << " (type: " << component.type << ", value: " << component.value << ")" << std::endl;  
        }  
    }

    // 打印节点  
    std::cout << "node:" << std::endl;  
    for (const auto &node : nodes) {  
        std::cout << node << std::endl;  
    }  

    // 打印边  
    std::cout << "\n edge:" << std::endl;  
    for (const auto &edge : edges) {  
        if (edge.type.empty()) {  
            std::cout << edge.component << ": " << edge.from << " -> " << edge.to   
                      << " (value: " << edge.value << ")" << std::endl;  
        } else {  
            std::cout << edge.component << ": " << edge.from << " -> " << edge.to   
                      << " (type: " << edge.type << ", value: " << edge.value << ")" << std::endl;  
        }  
    } 


    }
  return 1;
}
// cmake .. -Denable_testcases=ON -Denable_coverage=ON -DCMAKE_BUILD_TYPE=Debug -G "MSYS Makefiles"
// make