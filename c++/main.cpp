/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "hal.h"
#include "halCLParserInstance.h"

#include "psl_io.cpp"
#include "psl_merger.cpp"

static hal::CLParserPtr initParser() {
  hal::CLParserPtr optionsParser = hal::hdf5CLParserInstance(true);
  optionsParser->addArgument("inPslPath", "input psl file from liftover");
  optionsParser->addArgument("outPslPath", "output psl file ffor synteny blocks");
  //optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addArgument("minBlockSize", 
                           "lower bound on synteny block length. default 5000bp");
  optionsParser->addArgument("maxAnchorDistance", 
                           "upper bound on distance for syntenic psl blocks. default 5000bp");
  return optionsParser;
    
}


int main(int argc, char *argv[]) {
    //auto path = "/Users/admin/projects/breaks/data/AcinonyxJubatus.FelisCatus.19.psl";
    
//    auto min_block_size = 5000;
//    auto max_anchor_distance = 5000;
    hal::CLParserPtr optionsParser = initParser();
    std::string inPslPath;
    std::string outPslPath;
    int minBlockSize;
    int maxAnchorDistance;
    try {
        optionsParser->parseOptions(argc, argv);
        inPslPath = optionsParser->getArgument<std::string>("inPslPath");
        outPslPath = optionsParser->getArgument<std::string>("outPslPath");
        minBlockSize = optionsParser->getArgument<int>("minBlockSize");
        maxAnchorDistance = optionsParser->getArgument<int>("maxAnchorDistance");
        optionsParser->setDescription("convert psl alignments into synteny blocks");
    }
    catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        optionsParser->printUsage(std::cerr);
        exit(1);
    }
    auto blocks = psl_io::get_blocks_set(inPslPath);
    auto merged_blocks = dag_merge(blocks, minBlockSize, maxAnchorDistance);
    psl_io::write_psl(merged_blocks, outPslPath);
//    for (auto path : merged_blocks){
//        auto psl = psl_io::construct_psl(path);
//        std::cout << psl << std::endl;
//    } 
    return 0;
}