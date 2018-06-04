#include <string>
#include <numeric>
#include <sstream>
#include <iterator>

class PslBlock {
    //Block of a PSL
    public:
    int qStart;
    int qEnd;
    int tStart;
    int tEnd;
    int size;
    int tSize;
    int qSize;
    std::string qName;
    std::string tName;
    std::string tSeq;
    std::string qSeq;
    std::string strand;
    PslBlock(int qStart, int tStart, int size, std::string strand,
             const std::string& qName, const std::string& tName,
             int qSize, int tSize,
             const std::string& qSeq = "",
             const std::string& tSeq = ""){
        this->qStart = qStart;
        this->qSize = qSize;
        this->qEnd = qStart + size;
        this->tStart = tStart;
        this->tEnd = tStart + size;
        this->tSize = tSize;
        this->strand = strand;
        this->size = size;
        this->qName = qName;
        this->tName = tName;
        this->tSeq = tSeq;
        this->qSeq = qSeq;
    }
};


class Psl {
    public: 
    int match;
    int misMatch;
    int repMatch;
    int nCount;
    int qNumInsert;
    int qBaseInsert;
    int tNumInsert;
    int tBaseInsert;
    std::string strand;
    std::string qName;
    int qSize;
    int qStart;
    int qEnd;
    std::string tName;
    int tSize;
    int tStart;
    int tEnd;
    int blockCount;
    bool haveSeqs;
    std::vector<PslBlock> blocks;

    private:
    std::vector<std::string> split(const std::string &s, char delim) {
        std::stringstream ss(s);
        std::string item;
        std::vector<std::string> elems;
        while (std::getline(ss, item, delim)) {
            elems.push_back(std::move(item));
        }
        return elems;
    }
    
    std::vector<int> intArraySplit(const std::string& commaStr){
        //parser for comma-separated string list into a list of ints
        std::vector<int> ints;
        for (auto s : split(commaStr,','))
            ints.push_back(stoi(s));
        return ints;
    }

    void parseBlocks(const std::string& blockSizesStr,
                    const std::string& qStartsStr,
                    const std::string& tStartsStr,
                    int qSize, int tSize,
                    const std::string& qSeqsStr, 
                    const std::string& tSeqsStr) {
        auto blockSizes = intArraySplit(blockSizesStr);
        auto qStarts = intArraySplit(qStartsStr);
        auto tStarts = intArraySplit(tStartsStr);
        auto haveSeqs = (not qSeqsStr.empty());
        std::vector<std::string> qSeqs;
        std::vector<std::string> tSeqs;
        if (haveSeqs){
            qSeqs = split(qSeqsStr,',');
            tSeqs = split(tSeqsStr,',');
        }
        /*
        std::cout << "blockCount " << blockCount << std::endl;
        for (auto b : blockSizes) {
            std::cout << b << " ";
        }
        std::cout << std::endl;
        */
        for (int i = 0; i < this->blockCount; ++i){
            this->blocks.push_back(PslBlock( qStarts[i], tStarts[i], blockSizes[i], strand,
                                             qName, tName, qSize, tSize, 
                                         (haveSeqs ? qSeqs[i] : ""),
                                         (haveSeqs ? tSeqs[i] : "")));
        }
    }
    
    void parse(std::vector<std::string> row) {
        match = stoi(row[0]);
        misMatch = stoi(row[1]);
        repMatch = stoi(row[2]);
        nCount = stoi(row[3]);
        qNumInsert = stoi(row[4]);
        qBaseInsert = stoi(row[5]);
        tNumInsert = stoi(row[6]);
        tBaseInsert = stoi(row[7]);
        strand = row[8];
        qName = row[9];
        qSize = stoi(row[10]);
        qStart = stoi(row[11]);
        qEnd = stoi(row[12]);
        tName = row[13];
        tSize = stoi(row[14]);
        tStart = stoi(row[15]);
        tEnd = stoi(row[16]);
        blockCount = stoi(row[17]);
        haveSeqs = row.size() > 21;
        parseBlocks(row[18], row[19], row[20], qSize, tSize,
                           (haveSeqs ? row[21] : ""),
                           (haveSeqs ? row[22] : ""));
    }

    //template <class T>
    std::string vector_join(std::vector<std::string> v, std::string sep) const {
        std::string s = "";
        for (int i = 0; i < v.size(); ++i) {
            s += v[i];
            s += (i < v.size()-1 ? sep : "");
        }
        return s;
    }
    public:
    Psl(const std::string& line) {
        auto v = split(line,'\t');
        if (not (v.size() == 1) or line[0] == '#'){
            parse(v); 
        }
    }

    Psl() {}

    std::vector<PslBlock> get_blocks() {
        return blocks;
    }

    friend std::ostream& operator<<(std::ostream &strm, const Psl &a) ;

};

 std::ostream& operator<<(std::ostream &strm, const Psl &a) {    

        std::vector<std::string> v;
        v.push_back(std::to_string(a.match));
        v.push_back(std::to_string(a.misMatch));
        v.push_back(std::to_string(a.repMatch));
        v.push_back(std::to_string(a.nCount));
        v.push_back(std::to_string(a.qNumInsert));
        v.push_back(std::to_string(a.qBaseInsert));
        v.push_back(std::to_string(a.tNumInsert));
        v.push_back(std::to_string(a.tBaseInsert));
        v.push_back(a.strand);
        v.push_back(a.qName);
        v.push_back(std::to_string(a.qSize));
        v.push_back(std::to_string(a.qStart));
        v.push_back(std::to_string(a.qEnd));
        v.push_back(a.tName);
        v.push_back(std::to_string(a.tSize));
        v.push_back(std::to_string(a.tStart));
        v.push_back(std::to_string(a.tEnd));
        v.push_back(std::to_string(a.blockCount));
        std::vector<std::string> ints;
        std::vector<std::string> qstarts;
        std::vector<std::string> tstarts;
        for (auto b: a.blocks){
            ints.push_back(std::to_string(b.size));
            qstarts.push_back(std::to_string(b.qStart));
            tstarts.push_back(std::to_string(b.tStart));
        }
        v.push_back(a.vector_join(ints, ",")+",");
        v.push_back(a.vector_join(qstarts, ",")+",");
        v.push_back(a.vector_join(tstarts, ",")+",");
        //TODO add qseqs and tseqs!
        return strm << a.vector_join(v, "\t");
    }

