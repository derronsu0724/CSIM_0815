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
};  

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
        netlist << ".model NMOS NMOS (VTO=0.7)\n"; // NMOS模型定义  
        netlist << "M1 3 1 0 0 NMOS W=10u L=1u\n"; // 晶体管  
        netlist << ".tran 0.1ms 10ms\n";       // 转换分析指令  
        netlist << ".end\n";                   // 网表结束  

        // 关闭文件  
        netlist.close();  

        std::cout << "finish" << std::endl;  
    }
    else if  (temp1.at(1) == "4") {
      std::ifstream netlist("circuit.sp");  
      std::string line;  

      std::set<std::string> nodes;  // 用于存储节点  
      std::vector<Edge> edges;       // 用于存储边  

      if (!netlist) {  
          std::cerr << "无法打开文件！" << std::endl;  
          return 1;  
      }  

      while (std::getline(netlist, line)) {  
          // 跳过注释行  
          if (line.empty() || line[0] == '*') continue;  

          std::istringstream iss(line);  
          std::string component;  
          iss >> component; // 读取元件名称（例如 V1, R1）  

          std::string node1, node2;  
          // 解析节点  
          if (component[0] == 'V' || component[0] == 'R' || component[0] == 'C') {  
              iss >> node1 >> node2;  
          } else if (component[0] == 'M') { // 处理晶体管  
              iss >> node1 >> node2; // 读取输入两个节点  
              std::string otherNodes; // 其余的节点（如果有的话）  
              iss >> otherNodes; // 一般会有多个其它节点，且格式为: node3 node4...  
              // 你可以处理其它节点，如果需要的话  
          }  

          // 存储节点  
          nodes.insert(node1);  
          nodes.insert(node2);  

          // 存储边  
          edges.push_back({node1, node2, component});  
      }  

      netlist.close();  

      // 打印节点  
      std::cout << "node:" << std::endl;  
      for (const auto &node : nodes) {  
          std::cout << node << std::endl;  
      }  

      // 打印边  
      std::cout << "\n edge:" << std::endl;  
      for (const auto &edge : edges) {  
          std::cout << edge.component << ": " << edge.from << " -> " << edge.to << std::endl;  
      }
    }
  return 1;
}
