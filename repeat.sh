set -e
SCRIPT_PATH=$(realpath "$(dirname -- "${BASH_SOURCE[0]}")")
cd ${SCRIPT_PATH}/build/cpp/app/
rm -f ./*for_12_beams.txt
for run in {1..10}; do
    cp ${SCRIPT_PATH}/repetitions/${run}_6_5e-04/*moves_for_12_beams.txt .
    ./sc-runner 6 0.0005
    mkdir -p ${run}_6_5e-04
    mv ./field_* ./${run}_6_5e-04
    mv ./*for_12_beams.txt ./${run}_6_5e-04
done