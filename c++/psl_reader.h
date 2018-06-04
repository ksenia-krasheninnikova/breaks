#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "psl.h"


namespace psl_io {
std::vector<std::string> split(const std::string &s, char delim) {
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> elems;
  while (std::getline(ss, item, delim)) {
    elems.push_back(std::move(item)); 
  }
  return elems;
}


std::vector<PslBlock> get_blocks_set(const std::string psl) {
    std::ifstream input(psl);
    std::vector<PslBlock> blocks;
    for (std::string line; std::getline(input, line);) {
        auto cur_b = Psl(line).get_blocks();
        blocks.insert(blocks.end(), cur_b.begin(), cur_b.end());
    }
    return blocks;
 }

std::vector<int> get_qInserts(const std::vector<PslBlock>& blocks) {
    std::vector<int> result;
    for (int i=0; i < blocks.size()-1; ++i) {
        auto dif = blocks[i+1].qStart - blocks[i].qEnd;
        if (dif > 0) {
            //std::cout << dif << std::endl;
            result.push_back(dif);
        }
    }
    return result;
}

std::vector<int> get_tInserts(const std::vector<PslBlock>& blocks) {
    std::vector<int> result;
    for (int i=0; i < blocks.size()-1; ++i) {
        auto dif = blocks[i+1].tStart - blocks[i].tEnd;
        if (dif > 0) {
            result.push_back(dif);
        }
    }
    return result;
}

 Psl construct_psl(std::vector<PslBlock> blocks) {
    auto psl = Psl();
    psl.match = std::accumulate(blocks.begin(), blocks.end(), 0, [] (int total, PslBlock item) { return total + item.qEnd - item.qStart; });
    psl.misMatch = 0;
    psl.repMatch = 0;
    psl.nCount = 0;
    auto qInserts = get_qInserts(blocks);
    psl.qNumInsert = qInserts.size();
    psl.qBaseInsert = std::accumulate(qInserts.begin(), qInserts.end(), 0);
    auto tInserts = get_tInserts(blocks);
    psl.tNumInsert = tInserts.size();
    psl.tBaseInsert = std::accumulate(tInserts.begin(), tInserts.end(), 0);
    psl.qName = blocks[0].qName;
    psl.qSize = blocks[0].qSize;
    psl.qStart = blocks[0].qStart;
    psl.qEnd = blocks.back().qEnd;
    psl.tName = blocks[0].tName;
    psl.tSize = blocks[0].tSize;
    psl.strand = blocks[0].strand; 
    psl.tStart = blocks.back().tStart;
    psl.tEnd = blocks.back().tEnd;
    psl.blockCount = blocks.size();
    if (psl.strand == "++") {
        psl.tStart = blocks[0].tStart;
        psl.tEnd = blocks.back().tEnd;
    }
    else if (psl.strand == "+-") {
        psl.tEnd = psl.tSize - blocks[0].tStart;
        psl.tStart = psl.tSize - blocks.back().tEnd; 
    }
    psl.blocks = blocks;
    return psl;
 }
 }


