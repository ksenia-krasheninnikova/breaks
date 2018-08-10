synteny/ contains code for convertions of pairwise psl alignments into synteny blocks

gene_breaks/ contains code for convertions of single nucleotide markers into breakpoints

c++/ contains code for synteny implementaion as part of hal

Example usage:

Example:

Set ```rootPath``` in Makefile
```
cd c++
make
./hal2synteny /public/groups/cgl/felidae/cactus.hal out.psl --maxAnchorDistance 5000 --minBlockSize 5000 --queryGenome AcinonyxJubatus --targetGenome FelisCatus
```

Search for breakpoints in query.
```bash
${REPO_PATH}/order_query_by_target.py --reference ${REPO_PATH}/data/hits/A.bed --query ${REPO_PATH}/data/hits/B.bed | sort | uniq > ${REPO_PATH}/data/breaks/case1.breaks
```
Search for breakpoints in query related only to the speciefied seqs (chromosomes/scaffolds) in reference.
```bash
${REPO_PATH}/order_query_by_target.py --reference ${REPO_PATH}/data/hits/A.bed --query ${REPO_PATH}/data/hits/B.bed --seqs ${REPO_PATH}/data/seqs.txt | sort | uniq > ${REPO_PATH}/data/breaks/case2.breaks
```
In case of several breakpoints files sum up those breakpoints that are shared by different files:
```bash
${REPO_PATH}/sumup_breaks.py --breaks ${REPO_PATH}/data/breaks/ --hits ${REPO_PATH}/data/hits --ref A > ${REPO_PATH}/data/case_sumup.txt
```

Testing:
```bash
${REPO_PATH}/test.py
```
