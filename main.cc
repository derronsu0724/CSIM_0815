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
#include <spdlog/spdlog.h> // NOLINT
#include <spdlog/sinks/rotating_file_sink.h>  // NOLINT
#include "spdlog/sinks/stdout_color_sinks.h"  // NOLINT
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
        for (const auto &node : nodes) {  
            std::cout << node << ",";
        }
        std::cout << std::endl;*/
        // 打印边
        // std::cout << "\n edge:" << std::endl;
        for (const auto &edge : edges) {
            if (edge.type.empty()) {
                // std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (value: " << edge.value << ")" << std::endl;  
            } else if (edge.type == "DC") {            
                // std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_VDC);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "V", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::stod(edge.value)));
            } else if (edge.type == "resistor") {            
                // std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_R);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "R", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::stod(edge.value)));
            } else if (edge.type == "capacitor") {            
                // std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_CAP);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "C", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::stod(edge.value)));
            } else if (edge.type == "inductor") {            
                // std::cout << edge.component << ": " << edge.from << " -> " << edge.to << " (type: " << edge.type << ", value: " << edge.value << ")" << std::endl;
                ret = circuit->netlist()->addComponent(edge.component.c_str(), e_L);
                ret = circuit->netlist()->configComponent(edge.component.c_str(), "L", csimModel::Variant(csimModel::Variant::VariantDouble).setDouble(std::stod(edge.value)));
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

static int ACLinearCircuit(std::set<std::string> nodes,std::vector<Edge> edges,const char *fspace, double fstart, double fstop)
{
    std::cout   <<__LINE__ <<"ACLinearCircuit\n";
    int ret = 0;
    return ret;
}
void testSpdLog() {
	// 文件日志定义，设定日志最大100M，且最多保留10个
	auto fileLogger = spdlog::rotating_logger_mt("fileLogger", "log.log", 1024 * 1024 * 100, 10);
	//大于等于该等级的将被输出
	fileLogger->set_level(spdlog::level::trace);
 
	//定义输出格式
	spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l]  [%s->%!->%#]  [%v]");  //log类型 文件->函数名->行数   输出内容
	spdlog::set_default_logger(fileLogger);
 
	for (int i = 0; i < 2; ++i) {
		SPDLOG_LOGGER_TRACE(fileLogger, "trace输出：{}-{}", i, "测试trace");
		SPDLOG_LOGGER_DEBUG(fileLogger, "debug输出：{}-{}", i, "测试debug");
		SPDLOG_LOGGER_INFO(fileLogger, "info输出：{}-{}", i, "测试info");
		SPDLOG_LOGGER_WARN(fileLogger, "warn输出：{}-{}", i, "测试warn");
		SPDLOG_LOGGER_ERROR(fileLogger, "Some error message");
		SPDLOG_LOGGER_CRITICAL(fileLogger, "Some critical message");
	}
}

void testMultiLog() {
	//文件sink
	auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_st>("logs.log", 1024 * 1024 * 100, 10);
	file_sink->set_level(spdlog::level::debug);
	file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l]  [%s->%!->%#]  [%v]");
 
	//控制台sink
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
	console_sink->set_level(spdlog::level::debug);
	console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l]  [%s->%!->%#]  [%v]");
 
	std::vector<spdlog::sink_ptr> sinks;
	sinks.push_back(console_sink);
	sinks.push_back(file_sink);
 
	auto multiLogger = std::make_shared<spdlog::logger>("multiSink", begin(sinks), end(sinks));
	multiLogger->set_level(spdlog::level::debug);
	spdlog::set_default_logger(multiLogger);
 
	for (int i = 0; i < 2; ++i) {
		SPDLOG_LOGGER_TRACE(multiLogger, "trace输出：{}-{}", i, "测试trace");
		SPDLOG_LOGGER_DEBUG(multiLogger, "debug输出：{}-{}", i, "测试debug");
		SPDLOG_LOGGER_INFO(multiLogger, "info输出：{}-{}", i, "测试info");
		SPDLOG_LOGGER_WARN(multiLogger, "warn输出：{}-{}", i, "测试warn");
		SPDLOG_LOGGER_ERROR(multiLogger, "Some error message");
		SPDLOG_LOGGER_CRITICAL(multiLogger, "Some critical message");
	}
}



int main(int argc, char *argv[]) {
    // 创建控制台日志器  
    //auto console = spdlog::stdout_color_mt("console");
    // 创建一个文件日志器，指定日志文件的路径  
    // auto file_logger = spdlog::basic_logger_mt("file_logger", "logs/log.txt"); 
    //console->set_pattern("%^[%H:%M:%S] %v%$");  // 设置日志格式  
    //console->info("Hello, {}!", "world");  // 输出信息级日志  
    //console->warn("Warning message");       // 输出警告级日志  
    //console->error("Error occurred: {}", 404); // 输出错误级日志  
    // 使用默认日志器  
    auto console = spdlog::stdout_color_mt("console");    
    //auto err_logger = spdlog::stderr_color_mt("stderr");    
    spdlog::get("console")->info("loggers can be retrieved from a global registry using the spdlog::get(logger_name)");
    
	testMultiLog();
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
            std::string value; // 存储元件的值
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
        }*/

        OPCircuit_helper(nodes,edges, probe_names);
    }
  return 1;
}
// cmake .. -Denable_testcases=ON -Denable_coverage=ON -DCMAKE_BUILD_TYPE=Debug -G "MSYS Makefiles"
// make