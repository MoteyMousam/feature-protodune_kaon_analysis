source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source /storage/epp2/phrwdg/Dune/protoduneana/localProducts_larsoft_v09_12_00_e19_prof/setup
setup protoduneana v09_12_00 -q e19:prof
mrbsetenv
mrbslp
export FW_SEARCH_PATH=/storage/epp2/phrwdg/Dune/protoduneana:$FW_SEARCH_PATH
export FHICL_FILE_PATH=/storage/epp2/phrwdg/Dune/protoduneana:$FHICL_FILE_PATH
