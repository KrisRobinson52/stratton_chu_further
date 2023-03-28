set -e
SCRIPT_PATH=$(realpath "$(dirname -- "${BASH_SOURCE[0]}")")
cd ${SCRIPT_PATH}/build/cpp/app/
rm -f ./*for_12_beams.txt
for run in {1..10}; do
    ./sc-runner 6 0.0005
    mkdir -p ${run}_6_5e-04
    mv ./field_* ./${run}_6_5e-04
    mv ./*for_12_beams.txt ./${run}_6_5e-04
done