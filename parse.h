#ifndef CSIM_PARSE_H_
#define CSIM_PARSE_H_
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <map>
#include <variant>
namespace spice
{
    // 定义一个结构体来存储时间-电压对
    struct TimeVoltage {
        double time;
        double voltage;
    };
    struct Edge {
        std::string from;
        std::string to;
        std::string component;
        std::string type;
        //std::unordered_map<std::string, double> value; // 以名称为Key的参数值,电阻电容电感值
        std::unordered_map<std::string, std::variant<double,std::vector<double>, std::vector<TimeVoltage>>> value;
    };


    struct Subcircuit {
        std::string name; //  子电路名  
        std::vector<std::string> nodes;
        std::vector<std::string> external_nodes;
        std::vector<Edge> components; //  子电路内部元件  
    };
    double convert_to_hz(const std::string& value_str);
    std::vector<std::string> split_space(const std::string& text);
    double extractValues(const std::string& input, const std::string& keyword);
    std::vector<double> parseSIN(const std::string &str);
    std::vector<TimeVoltage> parsePWL(const std::string& str);
    int findNodePosition(const std::string& node1, const std::vector<std::string>& a1);
    std::map<std::string, std::vector<std::pair<spice::Edge, std::string>>> findEdgesForNodes(const std::vector<spice::Edge>& edges);
    void printEdgePairs(const std::map<std::string, std::vector<std::pair<spice::Edge, std::string>>>& nodeEdges);
    std::vector<std::pair<std::string, int>> findEdgesBetweenProbes(std::vector<spice::Edge> edges,std::vector<std::string> probe_names);
    void fun1(std::string& line,std::vector<std::string>& node,std::string& type, std::unordered_map<std::string, std::variant<double,std::vector<double>, std::vector<TimeVoltage>>>& data1);

    class parse_
    {
        public:
            parse_(std::string& s);
            ~parse_();
            std::vector<Subcircuit> m_subcircuits; // 存储所有子电路
            std::set<std::string> m_models;
            std::set<std::string> m_nodes;
            std::vector<Edge> m_edges;
            std::vector<std::string> m_probe_names;
            std::unordered_multimap<std::string, std::vector<double>>  m_analysis;            

    };
    /*
    class Graph {
        private:
            std::vector<std::vector<int>> adjacencyList;
            std::vector<bool> visited;
        public:
            Graph(int vertices);
            void addEdge(int u, int v);
            void DFS(int vertex);
    };*/
} // namespace spice
#endif // CSIM_PARSE_H_