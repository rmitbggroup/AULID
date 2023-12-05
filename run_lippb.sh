#!/usr/bin/env bash

BASE_PATH=/Users/lanhai/Program/learned_idx_disk_new/ALEX

DATASETPATH=/Users/lanhai/Program/data
PROGRAMS=(4k_disk)
DATASETS=(ycsb fb osm)

for FILE in ${PROGRAMS[@]}; do
  for DATASET in ${DATASETS[@]}; do
    HAS_SIZE=1
    TOTAL_COUNT=200000000
    W_TOTAL_COUNT=20000000

    if [ "$DATASET" = "ycsb" ];then
      HAS_SIZE=0
    fi

    if [ "$DATASET" = "osm800" ];then
      TOTAL_COUNT=800000000
      W_TOTAL_COUNT=200000000
    fi
    LOG1=${FILE}_${DATASET}.log1
    LOG2=${FILE}_${DATASET}.log2
    echo "Running ${FILE} on ${DATASET}"
    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=bsr.index --total_count=${TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=bulk_search_range --case_id=1 2>&1 >> ${BASE_PATH}/${LOG1}
    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=insert.index --total_count=${W_TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=insert 2>&1 >> ${BASE_PATH}/${LOG1}
    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=m1.index --total_count=${W_TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=mix_workload --case_id=1 2>&1 >> ${BASE_PATH}/${LOG1}
    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=m2.index --total_count=${W_TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=mix_workload --case_id=2 2>&1 >> ${BASE_PATH}/${LOG1}
    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=m3.index --total_count=${W_TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=mix_workload --case_id=3 2>&1 >> ${BASE_PATH}/${LOG1}

#    echo "Running ${FILE} on ${DATASET} again"
#    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=bsr.index --total_count=${TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=bulk_search_range --case_id=1 2>&1 >> ${BASE_PATH}/${LOG2}
#    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=insert.index --total_count=${TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=insert 2>&1 >> ${BASE_PATH}/${LOG2}
#    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=m1.index --total_count=${TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=mix_workload --case_id=1 2>&1 >> ${BASE_PATH}/${LOG2}
#    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=m2.index --total_count=${TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=mix_workload --case_id=2 2>&1 >> ${BASE_PATH}/${LOG2}
#    ${BASE_PATH}/${FILE} --keys_file=${DATASETPATH}/${DATASET} --index_file=m3.index --total_count=${TOTAL_COUNT} --has_size=${HAS_SIZE} --op_type=mix_workload --case_id=3 2>&1 >> ${BASE_PATH}/${LOG2}
  done
done