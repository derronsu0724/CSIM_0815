
#include "parse.h"
#include <sstream> 
#include <unordered_map>
#include <string>
#include <vector>
#include <boost/regex.hpp>

#include <iostream>
#include <fstream>
namespace spice
{
    std::vector<std::string> split_space(const std::string& text) {
    boost::regex ws_re("\\s+");
    std::vector<std::string> vector_(
        boost::sregex_token_iterator(text.begin(), text.end(), ws_re, -1),
        boost::sregex_token_iterator());
    return vector_;
    }

    double extractValues(const std::string& input, const std::string& keyword) {  
        double value= 0;  
        std::istringstream stream(input);  
        std::string temp;  

        while (stream >> temp) {  
            // 检查当前单词是否在关键词列表中  
            if (temp == keyword) {  
                stream >> value; // 获取关键词后面的值  
            }  
        }  

        return value;  
    }
    int findNodePosition(const std::string& node1, const std::vector<std::string>& a1) {  
        for (size_t i = 0; i < a1.size(); ++i) {  
            if (a1[i] == node1) {  
                return static_cast<int>(i); // 找到时返回索引  
            }  
        }  
        return -1; // -1 表示未找到  
    } 

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
              //auto a1=split_space(line);
              //node1=a1[1];
              //node2=a1[2];   
              iss >> node1 >> node2;        
              type="voltage";
              data1["DC"]=extractValues(line, "DC");
              data1["AC"]=extractValues(line, "AC");
          }
          node.push_back(node1);
          node.push_back(node2);
    }

    parse_::parse_(std::string& file) {
        std::cout<<file<<"\n";
        std::ifstream netlist(file);
        std::string line;
        // std::cout   <<__LINE__ <<"\n";
        if (netlist.is_open()) {
            Subcircuit current_subcircuit;  
            bool in_subcircuit = false;
            std::vector<Subcircuit> subcircuits; // 存储所有子电路
            
            while (std::getline(netlist, line)) {
                // std::cout << line << std::endl; // 打印每行内容  
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
                } else if (component == ".probe") {                
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
                                    m_probe_names.push_back(variable); // 将找到的变量保存到 vector 中
                                    start_pos = line.find("v(", end_pos);
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                    else if (component == ".op") {
                        m_analysis.emplace("op", std::vector<double>{1});  
                    }
                    else if (component == ".dc") {
                        auto a1=split_space(line);
                    }
                    else if (component == ".ac") {
                        auto a1=split_space(line);
                        std::vector<double> f0;
                        for(auto i=2;i<a1.size();i++){
                            f0.push_back(std::stod(a1[i]));
                        }
                        m_analysis.emplace("ac", f0);  
                    }
                    else if (component[0] == 'X') { 
                        //std::cout   <<__LINE__ <<"\n";
                        auto temp_=split_space(line);
                        std::string lastElement = temp_.back();
                        std::vector<std::string> a1;
                        for(size_t ii=1;ii<temp_.size()-1;ii++)
                        {
                            a1.push_back(temp_[ii]);
                            m_nodes.insert(temp_[ii]);
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
                                    m_edges.push_back({temp_node1, temp_node2, component.component+"_"+subcircuit.name, component.type, component.value});
                                }
                                break;
                            }
                            
                        }
                    }
            else{
                fun1(line,node_component,type,value);
                m_edges.push_back({node_component[0], node_component[1], component, type, value});
                m_nodes.insert(node_component[0]);
                m_nodes.insert(node_component[1]);
            }

            }
            m_probe_names.push_back("0");
            netlist.close(); // 关闭文件  
        } else {  
            std::cerr << "Unable to open file for reading." << std::endl;  
        }
    }

    parse_::~parse_()
    {
    }

    /*
        Graph::Graph(int n) {
            adjacencyList.resize(n);
            visited.resize(n, false);
        }
    
        void Graph::addEdge(int u, int v) {
            adjacencyList[u].push_back(v);
            adjacencyList[v].push_back(u);
        }
    
        void Graph::DFS(int vertex) {
            visited[vertex] = true;
            std::cout << vertex << " ";
    
            for (int i = 0; i < adjacencyList[vertex].size(); i++) {
                int neighbor = adjacencyList[vertex][i];
                if (!visited[neighbor]) {
                    DFS(neighbor);
                }
            }
        }
    */


} // namespace spice