#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>

#include "psl_reader.h"


//assume a.start < b.start
bool are_syntenic(const PslBlock& a, const PslBlock& b){
    return a.qEnd <= b.qStart &&
            a.tEnd <= b.tStart &&
             a.tName == b.tName &&
                a.strand == b.strand;
}

bool is_not_overlapping_ordered_pair(const PslBlock& a, const PslBlock& b, const int threshold=5000) {
    //std::cout << a.qStart << " " << b.qStart << " " << are_syntenic(a, b) << " " << a.tStart << " " << b.tStart << " "
    //    << (0 <= b.qStart - a.qEnd < threshold) << " " << (0 <=  b.tStart - a.tEnd < threshold) << " " << b.tStart - a.tEnd << " " << threshold  << std::endl;
    return are_syntenic(a, b) &&
             0 <= b.qStart - a.qEnd && 
                b.qStart - a.qEnd < threshold && 
                 0 <= b.tStart - a.tEnd &&
                    b.tStart - a.tEnd < threshold;
}

std::vector<int> get_next(const int pos, const std::vector<PslBlock>& queryGroup,
                 const int maxAnchorDistance=5000) {
    std::vector<int> f;
    //for (auto l : queryGroup) std::cout << l.qStart << " ";
    //std::cout << std::endl;
    for (auto i = pos + 1; i < queryGroup.size(); ++i) {
        if (is_not_overlapping_ordered_pair(queryGroup[pos], queryGroup[i], maxAnchorDistance)) {
            if (f.empty()) f.push_back(i);
            else {
                if (not f.empty() && is_not_overlapping_ordered_pair(queryGroup[f[0]], queryGroup[i], maxAnchorDistance)) {
                    //!!!!!!!!
                    ///for (auto l : f) std::cout << l << " (" << queryGroup[l].qStart << "," << queryGroup[l].qEnd << "," 
                    //                        << queryGroup[l].tStart << "," <<  queryGroup[l].tEnd << "," << queryGroup[l].strand << ") ";
                    //std::cout << std::endl;
                    return f;
                }
                else {
                    f.push_back(i);
                }
            }
        }
    }
    
    return f;
}

/*
dag is a dict that for a given vertex stores all its possible next verties
hidden_vertices is a set of vertices that are already in paths
weight of the edge equals length
of the next block
weight of a vertice equals
estimated weight: w_j < w_i + w_e(ij) =>
updated w_j
also remember how we came here:
(prev_vertex, weight)
*/
std::map<int, std::pair<int, int> > weigh_dag(const std::vector<PslBlock>& group, 
                                             std::map<int, std::vector<int> >& dag, 
                                             const std::set<int>& hiddenVertices,
                                             const int maxAnchorDistance){
    std::map<int,std::pair<int,int> > weightedDag;
    for (int i = 0; i < group.size(); ++i) {
        if (hiddenVertices.count(i)) continue; 
        std::vector<int> nexts;
        if (not dag.count(i)) {
            nexts = get_next(i, group, maxAnchorDistance);
            /*
            for (auto n : nexts) {
            std::cout << n << " ";
            }
            std::cout << std::endl;
            */
            dag[i] = nexts;
        }
        else {
            nexts = dag[i];
       }
       /*
       if never visited this vertex then
       its weight equals to its size
       because otherwise we will never count its size
       */
       if (not weightedDag.count(i))
            weightedDag[i] = std::make_pair(-1, group[i].size);
       for (auto j: nexts){
            if (hiddenVertices.count(j)) continue;
            auto alternativeWeight = weightedDag[i].second + group[j].size;
            if (not weightedDag.count(j) or weightedDag[j].second < alternativeWeight)
                //w_i + weight of the next edge 
                weightedDag[j] = std::make_pair(i, alternativeWeight);
       }
    }
    //std::cout << "weighted dag\n";
    //for (auto it = weightedDag.cbegin(); it != weightedDag.cend(); ++it) {
    //    std::cout << it ->first << " " << it->second.first << " " << it->second.second << std::endl;
    //}
    return weightedDag; 
}

int get_maxed_vertex(const std::map<int, std::pair<int, int> >& weightedDag) {
    int max = -10;
    int maxPos = -1;
    for (auto it = weightedDag.cbegin(); it != weightedDag.cend(); ++it) {
        //std::cout << "it->second.second " << it->second.second << std::endl;
        if (it->second.second >= max) {
            max = it->second.second;
            maxPos = it->first;
        }
    }
    //std::cout << "max pos " << maxPos << std::endl;
    return maxPos;
}

std::vector<PslBlock> traceback(std::map<int, std::pair<int, int> >& weightedDag, 
                                std::set<int>& hiddenVertices, 
                                const std::vector<PslBlock>& group) {
    //get the heaviest path weight
    auto startVertex = get_maxed_vertex(weightedDag);
    //std::cout << 11 << " " << startVertex << std::endl;
    std::vector<int> path = {startVertex};
    auto prevVertex = weightedDag[startVertex].first;
    //std::cout << 12 << std::endl;
    while (prevVertex != -1) {
        path.push_back(prevVertex);
        prevVertex = weightedDag[prevVertex].first;
    }
    hiddenVertices.insert(path.begin(), path.end());
    std::vector<PslBlock> pslBlockPath;
    //std::cout << "path size " << path.size() << std::endl;
    //for (auto p: path) {
    //    std::cout << p << " " << group[p].qStart << " ";
    //}
    //std::cout << std::endl;
    //std::cout << 4 << std::endl;
    std::transform(path.begin(), path.end(), std::back_inserter(pslBlockPath), [group] (int x) -> PslBlock {return group[x];});
    //std::cout << 5 << std::endl;
    std::reverse(pslBlockPath.begin(),pslBlockPath.end());
    //std::cout << 6 << std::endl;
    return pslBlockPath;
}
struct {    
    bool operator()(PslBlock a, PslBlock b) const
    {
        return a.qStart < b.qStart && a.tStart < b.tStart;
    }
} qStartLess;

std::vector<std::vector<PslBlock> >  dag_merge(const std::vector<PslBlock>& blocks, const int minBlockBreath, const int maxAnchorDistance){
    std::map<std::string, std::vector<PslBlock> > blocksByQName;
    for (auto block : blocks) blocksByQName[block.qName].push_back(block);
    std::vector<std::vector<PslBlock> > paths;
    for (auto pairs : blocksByQName) {
        std::vector<PslBlock> group = pairs.second; 
        std::map<int, std::vector<int> > dag;
        std::set<int> hiddenVertices;
        std::sort(group.begin(), group.end(), qStartLess); 
        while (hiddenVertices.size() != group.size()){
            //std::cout << 10 << std::endl;
            auto weightedDag = weigh_dag(group, dag, hiddenVertices, maxAnchorDistance);
            //std::cout << "group size " << group.size() << " hidden vertices size " << hiddenVertices.size() << std::endl;
            //std::cout << 2 << std::endl;
            auto path = traceback(weightedDag, hiddenVertices, group);
            //std::cout << path.size() << std::endl;
            if (path.empty()) {
                break;
            }
            //std::cout << 7 << std::endl;
            auto qLen = path.back().qEnd - path[0].qStart;
            auto tLen = path.back().tEnd - path[0].tStart;
            //std::cout << 8 << std::endl;
            if (qLen >= minBlockBreath && tLen >= minBlockBreath)
                paths.push_back(path);
            //std::cout << 9 << std::endl;

        }
    }
    return paths;
}


int main() {
    auto path = "/Users/admin/projects/breaks/data/AcinonyxJubatus.FelisCatus.19.psl";
    //auto path = "/hive/groups/recon/projs/felidae_comp/synteny-play-ground/data/felidae/cheetah/AcinonyxJubatus.FelisCatus.19.psl";
    auto min_block_size = 5000;
    auto max_anchor_distance = 5000;
    /*
    auto blocks = psl_io::get_blocks_set(path);
    std::cout << blocks.size() << "\n";
    for (auto b: blocks) {
        std::cout << b.qStart << " " << b.qEnd << " " << b.tStart << " " << b.tEnd << " " << b.size << "\n";
    }*/
    auto blocks = psl_io::get_blocks_set(path);
    auto merged_blocks = dag_merge(blocks, min_block_size, max_anchor_distance);
    for (auto path : merged_blocks){
        /*if (path[0].qName == "scaffold155" && path[0].qStart == 2646) { 
            for (auto b: path) {
                std::cout << b.qStart << " " << b.qEnd << " " << b.tStart << " " << b.tEnd << " " << b.size << " " << b.strand << "\n";
            }
            std::cout << path.size() << std::endl;
        }*/
        //if (path[0].qName == "scaffold155" && path[0].qStart == 2646) 
        auto psl = psl_io::construct_psl(path);
        std::cout << psl << std::endl;
    } 
    return 0;
}
