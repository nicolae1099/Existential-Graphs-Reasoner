// Copyright 2019 Burcea Marian-Gabriel, Nitu Nicolae-Iulian
#ifndef _HOME_MARIN_WORK_SD_TEMA3_AEGRAPH_H_
#define _HOME_MARIN_WORK_SD_TEMA3_AEGRAPH_H_

#include <vector>
#include <string>

class AEGraph {
 public:
    explicit AEGraph(std::string representation);

    std::string repr() const;

    void sort();

    bool operator<(const AEGraph& other) const;
    bool operator==(const AEGraph& other) const;
    bool operator!=(const AEGraph& other) const;
    AEGraph operator[](const int index) const;

    bool contains(const AEGraph& other) const;
    bool contains(const std::string other) const;

    int num_subgraphs() const;
    int num_atoms() const;
    int size() const;

    std::vector<std::vector<int>> possible_double_cuts();
    AEGraph double_cut(std::vector<int> where) const;

    std::vector<std::vector<int>> possible_erasures();
    AEGraph erase(std::vector<int>) const;
    void search_atom(std::string atom, AEGraph graph,
        std::vector<std::vector<int>>& paths);
    void search_graph(AEGraph to_find, AEGraph graph,
        std::vector<std::vector<int>>& paths);
    void set_paths(AEGraph& a);
    void set_levels(AEGraph& a);
    void set_cuts(AEGraph& a);
    std::vector<std::vector<int>> possible_deiterations();
    AEGraph deiterate(std::vector<int> where) const;
    std::vector<std::vector<int>> get_paths_to(const std::string other) const;
    std::vector<std::vector<int>> get_paths_to(const AEGraph& other) const;

    std::vector<std::string> atoms;
    std::vector<AEGraph> subgraphs;
    std::vector<int> path;
    bool level;
    bool cut;

    friend std::ostream& operator<<(std::ostream &out, const AEGraph &g);

    bool is_SA;
};

#endif  // _HOME_MARIN_WORK_SD_TEMA3_AEGRAPH_H_
