#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <ctime> 

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
            dag[i] = nexts;
        }
        else {
            nexts = dag[i];
       }
       //std::cout << i;
       //for (auto n : nexts) {
       //std::cout << " " << n;
       //}
       //std::cout << std::endl;
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
    //    std::cout << it ->first << " " << group[it->first].qStart <<  " " << group[it->first].qEnd << " " << group[it->first].tStart 
    //        << " "<< it->second.first << " " << it->second.second << std::endl;
    //}
    return weightedDag; 
}

int get_maxed_vertex(const std::map<int, std::pair<int, int> >& weightedDag) {
    int max = weightedDag.cbegin()->second.second;
    int maxPos = weightedDag.cbegin()->first;
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
    //    std::cout << group[p].qStart << " " << group[p].qEnd << " " << group[p].strand << " " << group[p].tStart << " " << group[p].tEnd << std::endl;
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
        //return a.qStart < b.qStart && a.tStart < b.tStart;
        if (a.qStart < b.qStart )
            return true;
        else if (a.qStart == b.qStart) {
            return a.tStart <= b.tStart;
        }
        return false;
    }
} qStartLess;

/*struct {    
    bool operator()(PslBlock a, PslBlock b) const
    {
        //return a.qStart < b.qStart && a.tStart < b.tStart;
        return a.tStart < b.tStart;
    }
} tStartLess;
*/

std::vector<std::vector<PslBlock> >  dag_merge(const std::vector<PslBlock>& blocks, const int minBlockBreath, const int maxAnchorDistance){
    std::map<std::string, std::vector<PslBlock> > blocksByQName;
    for (auto block : blocks) blocksByQName[block.qName].push_back(block);
    //for (auto pairs : blocksByQName) {
    //    std::cout << pairs.first << " " << pairs.second.size() << std::endl;
    //}
    std::vector<std::vector<PslBlock> > paths;
    for (auto pairs : blocksByQName) {
        //if (pairs.first != "scaffold1")
        //    continue;
        std::vector<PslBlock> group = pairs.second; 
        std::map<int, std::vector<int> > dag;
        std::set<int> hiddenVertices;
        std::sort(group.begin(), group.end(), qStartLess); 
        //std::stable_sort(group.begin(), group.end(), tStartLess); 
        //for (auto g: group){
        //    std::cout << g.qStart << " " << g.qName << " " << g.tStart << " " << g.tName<< std::endl;
        //}
        //std::cout << std::endl;
        while (hiddenVertices.size() != group.size()){
            //std::cout << "size " << hiddenVertices.size() << std::endl;
            //auto t = clock();
            auto weightedDag = weigh_dag(group, dag, hiddenVertices, maxAnchorDistance);
            //std::cout << "get_next time " << float(clock() - t)/CLOCKS_PER_SEC << std::endl;
            //t = clock();
            auto path = traceback(weightedDag, hiddenVertices, group);
            //std::cout << "traceback time " << float(clock() - t)/CLOCKS_PER_SEC << std::endl;
            //std::cout << path.size() << std::endl;
            if (path.empty()) {
                break;
            }
            auto qLen = path.back().qEnd - path[0].qStart;
            auto tLen = path.back().tEnd - path[0].tStart;
               // std::cout << path.size() << " " << tLen << " " << qLen << " " << minBlockBreath;
            if (qLen >= minBlockBreath && tLen >= minBlockBreath){
                paths.push_back(path);
                //std::cout << " OK" << std::endl;
            }
            //else
            //    std::cout << std::endl;
            //std::cout << pairs.first << " group size " << group.size() << " hidden vertices size " << hiddenVertices.size() << std::endl;
        }
    }
    return paths;
}


int main() {
    //auto path = "/Users/admin/projects/breaks/data/AcinonyxJubatus.FelisCatus.19.psl";
    //auto path = "/hive/groups/recon/projs/felidae_comp/synteny-play-ground/data/felidae/cheetah/AcinonyxJubatus.FelisCatus.19.psl";
    auto path = "/mnt/storagesm/projects/12cats/ksenia/felidae_comp/synteny-play-ground/data/felidae/cheetah/AcinonyxJubatus.FelisCatus.19.psl";
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
        auto psl = psl_io::construct_psl(path);
        std::cout << psl << std::endl;
    } 
    return 0;
}
