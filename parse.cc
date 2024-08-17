
#include "parse.h"
#include <sstream> 
#include <unordered_map>
#include <string>
#include <vector>
#include <boost/regex.hpp>
#include <variant>
#include <iostream>
#include <fstream>
#include <stdexcept>
namespace spice
{
    double convert_to_hz(const std::string& value_str) {
        double value = std::stod(value_str.substr(0, value_str.size() - 1));
        std::string unit = value_str.substr(value_str.size() - 1);

        if (unit == "k") {
            return value * 1000;
        } else if (unit == "M") {
            return value * 1000000;
        } else if (unit == "G") {
            return value * 1000000000;
        } else {
            throw std::invalid_argument("Invalid unit");
        }
    }

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

    std::vector<double> parseSIN(const std::string &str) {
        std::vector<double> result;
        std::istringstream iss(str);
        std::string token;

        // 跳过"SIN"之前的部分
        while (std::getline(iss, token, ' ')) {
            if (token == "sin") {
                break;
            }
        }

        // 解析"SIN"后面的数字
        while (std::getline(iss, token, ' ')) {
            try {
                //std::cout << token << " ";
                result.push_back(std::stod(token));
            } catch (const std::invalid_argument& e) {
                // 如果解析失败，跳过这个token
                continue;
            }
        }
        // 输出解析结果
        /*std::cout << "SIN\n";
        for (const auto& value : result) {
             std::cout << value << " ";
         }
        std::cout << "\n";*/
        return result;
    }

    // 解析PWL数据的函数
    std::vector<TimeVoltage> parsePWL(const std::string& str) {
        std::vector<TimeVoltage> data;
        std::string::size_type pos = str.find("pwl");
        if (pos != std::string::npos) {
            // 找到pwl后面的数据部分
            std::string pwlData = str.substr(pos + 3);
            std::istringstream iss(pwlData);
            std::string token;
            while (std::getline(iss, token, ',')) {
                // 分割每个时间-电压对
                std::istringstream tokenStream(token);
                double time, voltage;
                if (std::getline(tokenStream, token, 'n')) {
                    time = std::stod(token);
                    if (std::getline(tokenStream, token, 'v')) {
                        voltage = std::stod(token);
                        // 将时间-电压对添加到结果中
                        data.push_back({time, voltage});
                    }
                }
            }
        }
        // 输出解析结果
        // for (const auto& pair : data) {
        //     std::cout << "Time: " << pair.time << "n, Voltage: " << pair.voltage << "v" << std::endl;
        // }
        return data;
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
        /*
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
        */
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
                    const spice::Edge& edge1 = edgesInfo[i].first;  
                    const spice::Edge& edge2 = edgesInfo[j].first;
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
        /*
        for (const auto& ii  : componentFromPair) {
                    std::cout << "Component: " << ii.first << std::endl;  
                    std::cout << "node: " << ii.second << std::endl;
        }*/
        return componentFromPair;
    }

    void fun1(std::string& line,std::vector<std::string>& node,std::string& type, std::unordered_map<std::string, std::variant<double,std::vector<double>, std::vector<TimeVoltage>>>& data1)
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
              //data1["pwl"]=parsePWL(line);
              data1["sin"]=parseSIN(line);
                // 输出解析结果
                std::cout << "SIN\n";
                for (const auto& value : std::get<std::vector<double>>(data1["sin"])) {
                    std::cout << value << " ";
                }
                std::cout << "\n";
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
           
            bool in_subcircuit = false;
            std::vector<Subcircuit> subcircuits; // 存储所有子电路
            
            while (std::getline(netlist, line)) {
                // std::cout << line << std::endl; // 打印每行内容  
                // 跳过注释行  
                if (line.empty() || line[0] == '*') continue;
                Subcircuit current_subcircuit;
                std::istringstream iss(line);
                std::string component;
                iss >> component; // 读取元件名称
                std::vector<std::string> node_component;
                std::string type;  // 存储类型  
                //std::vector<double> value; // 存储元件的值
                //std::unordered_map<std::string, double> value; // 存储元件的值
                std::unordered_map<std::string, std::variant<double,std::vector<double>, std::vector<TimeVoltage>>> value; // 存储元件的值
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
                    else if (component == ".tran") {
                        std::string temp1;
                        std::vector<double> t0;
                        while (iss >> temp1) {
                            //std::cout <<temp1<<"\n";
                            t0.push_back(std::stod(temp1));
                            std::cout <<t0[t0.size()-1]<<"\n";
                        }
                        m_analysis.emplace("tran", t0);  
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
                    // 输出解析结果
                    /*std::unordered_map<std::string, double> value2;
                    for (const auto& pair : value) {
                        std::cout << "Key: " << pair.first << std::endl;
                        if (std::holds_alternative<double>(pair.second)) {
                            value2[pair.first]    =std::get<double>(pair.second) ;
                            std::cout << "Value: " << std::get<double>(pair.second) << std::endl;
                        } else if (std::holds_alternative<std::vector<TimeVoltage>>(pair.second)) {
                            for (const auto& tv : std::get<std::vector<TimeVoltage>>(pair.second)) {
                                std::cout << "Time: " << tv.time << ", Voltage: " << tv.voltage << std::endl;
                            }
                        }
                    }
                    */
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