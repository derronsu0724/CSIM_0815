#ifndef CSIM_PARSE_H_
#define CSIM_PARSE_H_
#include <vector>
#include <set>
#include <string>
#include <unordered_map>
namespace spice
{
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
    std::vector<std::string> split_space(const std::string& text);
    double extractValues(const std::string& input, const std::string& keyword);
    int findNodePosition(const std::string& node1, const std::vector<std::string>& a1);
    void fun1(std::string& line,std::vector<std::string>& node,std::string& type, std::unordered_map<std::string, double>& data1);

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