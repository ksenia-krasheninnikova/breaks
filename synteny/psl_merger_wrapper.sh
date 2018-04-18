#!/bin/bash
source ~/.bashrc
set -beEu -o pipefail
/hive/groups/recon/projs/felidae_comp/synteny-play-ground/bin/psl_merger.py dag 100000 100000 "$1" "$2"
