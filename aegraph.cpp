// Copyright 2019 Burcea Marian-Gabriel, Nitu Nicolae-Iulian
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <cassert>
#include <queue>
#include <stack>
#include "/home/marin/work/sd/tema3/aegraph.h"

void print_vector(std::vector<int> v) {
    for (unsigned int i = 0; i < v.size(); i++) {
        std::cout << v[i] << ' ';
    }
    std::cout << "\n\n";
}

void print_matrix(std::vector<std::vector<int>> v) {
    for (unsigned int i = 0; i < v.size(); i++) {
        for (unsigned int j = 0; j < v[i].size(); j++) {
            std::cout << v[i][j] << ' ';
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; i++) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }

    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; i++) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; i++) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::possible_double_cuts() {
    std::vector<std::vector<int>> paths;
    std::queue<AEGraph> q;
    AEGraph a = *this;
    set_paths(a);
    bool pos = false;
    q.push(a);
    while (!q.empty()) {
        AEGraph temp = q.front();
        q.pop();
        if (temp.num_subgraphs() == 1 && temp.num_atoms() == 0 && pos == true) {
            paths.push_back(temp.path);
        }
        for (int i = 0; i < temp.num_subgraphs(); i++) {
            q.push(temp.subgraphs[i]);
        }
        pos = true;
    }
    return paths;
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
    AEGraph p = *this;
    AEGraph* copy = &p;
    for (unsigned int i = 0; i < where.size() - 1; i++) {
        AEGraph* tmp = &copy->subgraphs[where[i]];
        copy = tmp;
    }
    AEGraph tmp = copy->subgraphs[where[where.size() - 1]];
    AEGraph x = tmp;
    AEGraph x2 = x;
    tmp = x.subgraphs[0];
    x = tmp;
    std::vector<AEGraph> sub;
    for (unsigned int i = 0; i < x.subgraphs.size(); i++) {
        sub.push_back(x.subgraphs[i]);
    }
    std::vector<std::string> at;
    for (unsigned int i = 0; i < x.atoms.size(); i++) {
        at.push_back(x.atoms[i]);
    }
    for (unsigned int i = 0; i < sub.size(); i++) {
        copy->subgraphs.push_back(sub[i]);
    }
    for (unsigned int i = 0; i < at.size(); i++) {
        copy->atoms.push_back(at[i]);
    }
    std::vector<AEGraph>::iterator it =
    find(copy->subgraphs.begin(), copy->subgraphs.end(), x2);
    copy->subgraphs.erase(it);
    return p;
}

void AEGraph::set_levels(AEGraph& a) {
    AEGraph* copy = &a;
    std::queue<AEGraph*> q;
    q.push(copy);
    copy->level = 0;
    while (!q.empty()) {
        AEGraph* x = q.front();
        q.pop();
        for (unsigned int i = 0; i < x->subgraphs.size(); i++) {
            x->subgraphs[i].level = 1 - x->level;
            q.push(&x->subgraphs[i]);
        }
    }
}

void AEGraph::set_cuts(AEGraph& a) {
    AEGraph* copy = &a;
    std::queue<AEGraph*> q;
    q.push(copy);
    copy->cut = 0;
    while (!q.empty()) {
        AEGraph* x = q.front();
        q.pop();
        for (unsigned int i = 0; i < x->subgraphs.size(); i++) {
            if (x->subgraphs.size() == 1 && x->atoms.size() == 0) {
                x->subgraphs[i].cut = 0;
            } else {
                x->subgraphs[i].cut = 1;
            }
            q.push(&x->subgraphs[i]);
        }
    }
}

std::vector<std::vector<int>> AEGraph::possible_erasures() {
    std::vector<std::vector<int>> paths;
    AEGraph a = *this;
    set_paths(a);
    set_levels(a);
    set_cuts(a);
    std::queue<AEGraph> q;
    for (unsigned int i = 0; i < a.subgraphs.size(); i++) {
        paths.push_back(a.subgraphs[i].path);
    }
    for (unsigned int i = 0; i < a.atoms.size(); i++) {
        std::vector<int> v;
        v.push_back(a.subgraphs.size() + i);
        paths.push_back(v);
    }
    q.push(a);
    while (!q.empty()) {
        AEGraph x = q.front();
        q.pop();
        for (unsigned int i = 0; i < x.subgraphs.size(); i++) {
            if (x.subgraphs[i].level == 1 && x.subgraphs[i].cut == 1) {
                bool exist = false;
                for (unsigned int j = 0; j < paths.size(); j++) {
                    if (paths[j] == x.subgraphs[i].path) {
                        exist = true;
                        break;
                    }
                }
                if (exist == false) {
                    paths.push_back(x.subgraphs[i].path);
                }
            } else if (x.subgraphs[i].level == 0) {
                if (x.subgraphs[i].atoms.size() > 1) {
                    for (unsigned int t = 0; t <
                        x.subgraphs[i].atoms.size(); t++) {
                        std::vector<int> v;
                        for (unsigned int j = 0;
                            j < x.subgraphs[i].path.size(); j++) {
                            v.push_back(x.subgraphs[i].path[j]);
                        }
                        v.push_back(t);
                        paths.push_back(v);
                    }
                }
            }
            q.push(x.subgraphs[i]);
        }
    }
    return paths;
}


AEGraph AEGraph::erase(std::vector<int> where) const {
    AEGraph p = *this;
    AEGraph* copy = &p;
    for (unsigned int i = 0; i < where.size() - 1; i++) {
        AEGraph* tmp = &copy->subgraphs[where[i]];
        copy = tmp;
    }
    AEGraph tmp = *copy;
    unsigned int last_pos = where[where.size() - 1];
    if (last_pos < copy->subgraphs.size()) {
        copy->subgraphs.erase(copy->subgraphs.begin() + last_pos);
    } else {
        last_pos = last_pos - copy->subgraphs.size();
        copy->atoms.erase(copy->atoms.begin() + last_pos);
    }
    return p;
}

void AEGraph::set_paths(AEGraph& a) {
    std::queue<AEGraph*> q;
    AEGraph* copy = &a;
    q.push(copy);
    while (!q.empty()) {
        AEGraph* x = q.front();
        q.pop();
        for (unsigned int i = 0; i < x->subgraphs.size(); i++) {
            for (unsigned int j = 0; j < x->path.size(); j++) {
                x->subgraphs[i].path.push_back(x->path[j]);
            }
            x->subgraphs[i].path.push_back(i);
            q.push(&x->subgraphs[i]);
        }
    }
}

void AEGraph::search_atom(std::string atom, AEGraph graph,
    std::vector<std::vector<int>>& paths) {
    std::queue<AEGraph> q;
    q.push(graph);
    while (!q.empty()) {
        AEGraph x = q.front();
        q.pop();
        for (unsigned int i = 0; i < x.atoms.size(); i++) {
            if (x.atoms[i] == atom) {
                std::vector<int> v;
                for (unsigned int j = 0; j < x.path.size(); j++) {
                    v.push_back(x.path[j]);
                }
                v.push_back(x.subgraphs.size() + i);
                bool exist = false;
                for (unsigned int k = 0; k < paths.size(); k++) {
                    if (paths[k] == v) {
                        exist = true;
                        break;
                    }
                }
                if (exist == false) {
                    paths.push_back(v);
                }
            }
        }
        if (x.subgraphs.size() == 1 && x.subgraphs[0].atoms.size() == 1
            && x.subgraphs[0].subgraphs.size() == 0) {
            continue;
        }
        for (unsigned int i = 0; i < x.subgraphs.size(); i++) {
            q.push(x.subgraphs[i]);
        }
    }
}

void AEGraph::search_graph(AEGraph to_find, AEGraph graph,
    std::vector<std::vector<int>>& paths) {
    std::queue<AEGraph> q;
    q.push(graph);
    while (!q.empty()) {
        AEGraph x = q.front();
        q.pop();
        for (unsigned int i = 0; i < x.subgraphs.size(); i++) {
            if (x.subgraphs[i].repr() == to_find.repr()) {
                std::vector<int> v;
                bool exist = false;
                for (unsigned int k = 0; k < paths.size(); k++) {
                    if (paths[k] == x.subgraphs[i].path) {
                        exist = true;
                    }
                }
                if (exist == false) {
                    paths.push_back(x.subgraphs[i].path);
                }
            }
            q.push(x.subgraphs[i]);
        }
    }
}

std::vector<std::vector<int>> AEGraph::possible_deiterations() {
    std::vector<std::vector<int>> paths;
    AEGraph a = *this;
    set_paths(a);
    std::queue<AEGraph> q;
    q.push(a);
    while (!q.empty()) {
        AEGraph x = q.front();
        q.pop();
        for (unsigned int i = 0; i < x.subgraphs.size(); i++) {
            for (unsigned int j = 0; j < x.atoms.size(); j++) {
                search_atom(x.atoms[j], x.subgraphs[i], paths);
            }
            for (unsigned int j = 0; j < x.subgraphs.size(); j++) {
                search_graph(x.subgraphs[j], x.subgraphs[i], paths);
            }
            q.push(x.subgraphs[i]);
        }
    }
    return paths;
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
    AEGraph p = *this;
    AEGraph* copy = &p;
    for (unsigned int i = 0; i < where.size() - 1; i++) {
        AEGraph* tmp = &copy->subgraphs[where[i]];
        copy = tmp;
    }
    AEGraph tmp = *copy;
    unsigned int last_pos = where[where.size() - 1];
    if (last_pos < copy->subgraphs.size()) {
        copy->subgraphs.erase(copy->subgraphs.begin() + last_pos);
    } else {
        last_pos = last_pos - copy->subgraphs.size();
        copy->atoms.erase(copy->atoms.begin() + last_pos);
    }
    return p;
}
