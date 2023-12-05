#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <random>
#include <chrono>
#include <map>
#include "storage_management.h"
#include <unistd.h>
#include <string>
#include "util.h"
#include "utils.h"

#define KeyType uint64_t
#define ValueType uint64_t
using namespace std;

void test_blipp_bulk(int memory_type, char *index_name, char *key_path, int count, int has_size) {
    LIPPBTree<KeyType, ValueType> index;
    index.init(index_name, true, memory_type);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);

    ValueType *values = new ValueType[count];
    for (int i = 0; i < count; i++) {
        values[i] = keys[i] + 1;
    }
    cout << "start to build... " << endl;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    index.bulk_load_entry(keys, values, count);
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    std::cout << "bulk load time: " << bulk_lookup_time / 1e9 << std::endl;
    std::cout << "file size:" << index.report_file_size() << " bytes" << std::endl;
    delete[]keys;
    delete[]values;
    return;
}

void
test_blipp_lookup(int memory_type, char *index_name, char *key_path, int count, int has_size, int s_count, int case_id,
                  int step, int need_sleep) {
    LIPPBTree<KeyType, ValueType> index;
    index.init(index_name, false, memory_type);
//    fstream file;
//    file.open("search_keys.txt");
//    long a;
//    file >> a;
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);
//
//    ValueType *values = new ValueType[count];
//    for (int i = 0; i < count; i++) {
//        values[i] = keys[i] + 1;
//    }
//    cout << "start bulk... " << endl;
////    index.bulk_load_entry(keys, values, count);
//    cout << "end bulk... " << endl;
    bool found;
    int bc;
    double sc = 0;
    ValueType v;
    bc = 0;
    KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
//    KeyType sk = 13748238711574469023;
    delete[] keys;

    cout << "start to test... " << endl;
    if (need_sleep == 1) sleep(7);
    cout << "start to record... " << endl;
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    int nc = 0;
    double in_t = 0;
    int ic = 0;

    for (int j = 0; j < s_count; j++) {
        bc = 0;
        ic = 0;
        found = index.lippb_search_entry(search_keys[j], &v, &bc, &ic);
        sc += bc;
        in_t += ic;
        if (!found) {
            std::cout << "not" << j << std::endl;
            return;
            nc += 1;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    cout << 1e9 * s_count / batch_lookup_time << " ops" << endl;
    cout << sc / s_count << " block/lookup" << endl;
    cout << in_t / s_count << " in block/lookup" << endl;
    cout << "not found:" << nc << endl;
    delete[]search_keys;
}

void
test_blipp_scan(int memory_type, char *index_name, char *key_path, int count, int has_size, int s_count, int case_id,
                int r_size, int step) {
    LIPPBTree<KeyType, ValueType> index;
    index.init(index_name, true, memory_type);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);

    ValueType *values = new ValueType[count];
    for (int i = 0; i < count; i++) {
        values[i] = keys[i] + 1;
    }

    cout << "start bulk... " << endl;
    index.bulk_load_entry(keys, values, count);
    cout << "end bulk... " << endl;

    bool found;
    int bc;
    double sc = 0;
    KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
    delete[] keys;
    ValueType v;
    KeyType rs[r_size];
    cout << "start to test... " << endl;

    cout << "start to scan... " << endl;
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    int nc = 0;
    double in_t = 0;
    int ic = 0;
    for (int j = 0; j < s_count; j++) {
        bc = 0;
        ic = 0;
        index.lippb_scan_entry(search_keys[j], rs, &bc, r_size, &ic);
        sc += bc;
        in_t += ic;
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    cout << 1e9 * s_count / batch_lookup_time << " ops" << endl;
    cout << sc / s_count << " block/lookup" << endl;
    cout << in_t / s_count << " in block/lookup" << endl;
    delete[]search_keys;
}

int get_next(std::vector<int> &v, int _seed) {
    int n = v.size();

    srand(_seed);

    // Make sure the number is within
    // the index range
    int index = rand() % n;

    // Get random number from the vector
    int num = v[index];

    // Remove the number from the vector
    std::swap(v[index], v[n - 1]);
    v.pop_back();
    return num;
}

void
test_insert(int memory_type, char *index_name, char *key_path, int count, int has_size, int insert_count = 10000000) {
    LIPPBTree<KeyType, ValueType> index;
    index.init(index_name, true, memory_type);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    // To obtain a time-based seed
    unsigned seed = 0;

    // Shuffling our array
    shuffle(keys, keys + count,
            default_random_engine(seed));
    std::sort(keys, keys + count - insert_count);
    ValueType *values = new ValueType[count - insert_count];
    for (int i = 0; i < count - insert_count; i++) {
        values[i] = keys[i] + 1;
    }
    int bc = 0;
    int nc = 0;
    double sc = 0;
    bool found;
    ValueType _v;
    double t = 0;

    cout << "start to bulk... " << endl;
    index.bulk_load_entry(keys, values, count - insert_count);
    cout << "finish to bulk... " << endl;
    std::vector<int> v(insert_count);
    for (int i = 0; i < insert_count; i++)
        v[i] = i;
    std::vector<int> i_o(insert_count);
    for (int i = 0; i < insert_count; i++) {
        i_o[i] = get_next(v, insert_count - i);
    }
    int _start = count - insert_count;
//    long long *cc_array = new long long[insert_count];
    double total_level = 0;
    int ic = 0;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    for (int i = count - insert_count, j = 0; i < count; i++, j++) {
        bc = 0;
//        std::chrono::high_resolution_clock::time_point i_s = std::chrono::high_resolution_clock::now();

        index.insert_lippb_entry(keys[_start + i_o[j]], keys[_start + i_o[j]] + 1, &bc, &ic);
//        std::chrono::high_resolution_clock::time_point i_e = std::chrono::high_resolution_clock::now();
//        cc_array[j] = std::chrono::duration_cast<std::chrono::nanoseconds>(i_e - i_s).count();
        t += bc;
        total_level += ic;

    }
    NodeHeaderD _head;
    index.sys_metablock(false, _head);
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    cout << double(insert_count) * 1e9 / bulk_lookup_time << " ops" << endl;
    cout << t / insert_count << " block/lookup" << endl;
    cout << total_level / insert_count << " node/lookup" << endl;
    std::cout << "file size:" << index.report_file_size() << " bytes" << std::endl;

    index.print_mm();
    index.print_mb();
    index.print_block();
    delete[]keys;
    delete[]values;
    return;
}

void hybrid_test(int memory_type, char *index_name, char *key_path, int count, int has_size, int case_id,
                 int total_op = 10000000) {
    LIPPBTree<KeyType, ValueType> index;
    index.init(index_name, true, memory_type);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    unsigned seed = 0;

    // Shuffling our array
    shuffle(keys, keys + count,
            default_random_engine(seed));
    std::sort(keys, keys + count - total_op);
    ValueType *values = new ValueType[count - total_op];
    for (int i = 0; i < count - total_op; i++) {
        values[i] = keys[i] + 1;
    }
//    int total_op = count/2;
    int start = count - total_op;
    std::vector<int> v(total_op);
    for (int i = 0; i < total_op; i++)
        v[i] = i;
    std::vector<KeyType> i_o(total_op);
    for (int i = 0; i < total_op; i++) {
        i_o[i] = keys[start + get_next(v, total_op - i)];
    }
    for (int i = 0; i < total_op; i++) {
        keys[start + i] = i_o[i];
    }
    int s_count;
    int step_s;
    int step_w;
    if (case_id == 1) { // read_heavey
        s_count = 18 * total_op / 20;
        step_s = 18;
        step_w = 2;
    } else if (case_id == 2) { // write_heavey
        s_count = 2 * total_op / 20;
        step_s = 2;
        step_w = 18;
    } else { // balence
        s_count = total_op / 2;
        step_s = 10;
        step_w = 10;
    }
    KeyType *search_keys = new KeyType[s_count];
    std::mt19937_64 gen(19937);
    for (int i = 0; i < s_count;) {
        std::uniform_int_distribution<int> dis(0, start + (i / step_s + 1) * step_w);
        for (int _i = 0; _i < step_s; _i++) {
            search_keys[i] = keys[dis(gen)];
            i++;
        }
    }

    int w_count = total_op - s_count;
    int bc = 0;
    bool found;
    ValueType _v;
    double t = 0;
    int nc = 0;
    cout << "start to bulk... " << endl;
    index.bulk_load_entry(keys, values, count - total_op);
    cout << "start to hybrid test... " << endl;
    int _start = count - total_op;
    double total_ic = 0;
    int ic = 0;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    for (int i = _start, j = 0; i < (_start + w_count) && j < s_count;) {
        for (int _i = 0; _i < step_w; _i++) {

            bc = 0;
            bool x = index.insert_lippb_entry(keys[i], keys[i] + 1, &bc, &ic);

            t += bc;
            total_ic += ic;
            i += 1;
        }
        for (int _j = 0; _j < step_s; _j++) {
//            std::cout<< "s-" << j << "-" <<search_keys[j] <<std::endl;
            bc = 0;
            int ic = 0;
            found = index.lippb_search_entry(search_keys[j], &_v, &bc, &ic);
            t += bc;
            j += 1;
            if (!found) {
//                std::cout << "ni-" << i << std::endl;
//                return;
                nc += 1;
            }
        }
    }
    NodeHeaderD _head;
    index.sys_metablock(false, _head);
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    double scount = total_op;
    cout << 1e9 * scount / bulk_lookup_time << " ops" << endl;
    cout << t / scount << " block/lookup" << endl;
    cout << nc << endl;
    std::cout << "file size:" << index.report_file_size() << " bytes" << std::endl;

    delete[]keys;
    delete[]values;
    delete[]search_keys;
    return;
}

void test_bulk_search(int memory_type, char *index_name, char *key_path, int count, int has_size,
                      int s_count, int case_id, int r_size) {
    LIPPBTree<KeyType, ValueType> index;
    index.init(index_name, true, memory_type);
    std::ifstream fin(key_path, std::ios::binary);
    KeyType *keys = new KeyType[count];
    if (has_size == 1) {
        uint64_t size;
        fin.read((char *) (&size), sizeof(uint64_t));
    }
    fin.read((char *) (keys), sizeof(KeyType) * count);
    fin.close();
    std::sort(keys, keys + count);

    ValueType *values = new ValueType[count];
    for (int i = 0; i < count; i++) {
        values[i] = keys[i] + 1;
    }
    cout << "start to build... " << endl;
    std::chrono::high_resolution_clock::time_point bulk_start = std::chrono::high_resolution_clock::now();
    index.bulk_load_entry(keys, values, count);
    std::chrono::high_resolution_clock::time_point bulk_end = std::chrono::high_resolution_clock::now();
    long long bulk_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(bulk_end - bulk_start).count();
    std::cout << "bulk load time: " << bulk_lookup_time / 1e9 << std::endl;
    std::cout << "file size:" << index.report_file_size() << " bytes" << std::endl;
    std::cout << "\n\n\n\n" << std::endl;
    delete[] values;


    bool found;
    int bc;
    double sc = 0;
    ValueType v;
    bc = 0;
    KeyType *search_keys;
    if (case_id == 1) {
        search_keys = get_search_keys(keys, count - 2 * r_size, s_count);
    } else if (case_id == 2) {
        search_keys = get_search_keys_zipf(keys, count - 2 * r_size, s_count);
    } else {
        throw std::invalid_argument("not support this query case...");
    }
    delete[] keys;
//    long long *latency = new long long[s_count];
    std::chrono::high_resolution_clock::time_point lookups_start_time = std::chrono::high_resolution_clock::now();
    int nc = 0;
    double in_t = 0;
    int ic = 0;
    for (int j = 0; j < s_count; j++) {
        bc = 0;
        ic = 0;
//        std::chrono::high_resolution_clock::time_point i_lookups_start_time = std::chrono::high_resolution_clock::now();
        found = index.lippb_search_entry(search_keys[j], &v, &bc, &ic);
//        std::chrono::high_resolution_clock::time_point i_lookups_end_time = std::chrono::high_resolution_clock::now();
//        long long i_batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
//                i_lookups_end_time - i_lookups_start_time).count();
//        latency[j] = i_batch_lookup_time;
        sc += bc;
        in_t += ic;
        if (!found) {
            nc += 1;
        }
    }
    std::chrono::high_resolution_clock::time_point lookups_end_time = std::chrono::high_resolution_clock::now();
    long long batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    cout << 1e9 * s_count / batch_lookup_time << " ops" << endl;
    cout << sc / s_count << " block/lookup" << endl;
    cout << in_t / s_count << " in block/lookup" << endl;
    cout << "not found:" << nc << endl;
    std::cout << "\n\n\n\n" << std::endl;
//    std::ofstream outFile;
//    outFile.open("s_latency.txt");
//    for (int i = 0; i < s_count; i++) {
//        outFile << latency[i] << std::endl;
//    }
//    outFile.close();
//    return;
    lookups_start_time = std::chrono::high_resolution_clock::now();
    nc = 0;
    in_t = 0;

    ic = 0;
    sc = 0;
    KeyType rs[r_size];
    for (int j = 0; j < s_count; j++) {
        bc = 0;
        ic = 0;
        index.lippb_scan_entry(search_keys[j], rs, &bc, r_size, &ic);
        sc += bc;
        in_t += ic;
    }
    lookups_end_time = std::chrono::high_resolution_clock::now();
    batch_lookup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            lookups_end_time - lookups_start_time).count();
    cout << 1e9 * s_count / batch_lookup_time << " ops" << endl;
    cout << sc / s_count << " block/lookup" << endl;
    cout << in_t / s_count << " in block/lookup" << endl;
    delete[]search_keys;

    std::cout << "\n\n\n\n" << std::endl;
}

int main(int argc, char *argv[]) {
    auto flags = parse_flags(argc, argv);
    std::string key_file_path = get_required(flags, "keys_file");
    std::string op_type = get_required(flags, "op_type");
    std::string index_name = get_required(flags, "index_file");
    int count = stoi(get_required(flags, "total_count"));
    int search_count = stoi(get_with_default(flags, "search_count", "200000"));
    int has_size = stoi(get_required(flags, "has_size"));
    int case_id = stoi(get_with_default(flags, "case_id", "1"));
    int step = stoi(get_with_default(flags, "step", "999"));
    int r_size = stoi(get_with_default(flags, "r_size", "100"));
    int need_sleep = stoi(get_with_default(flags, "need_sleep", "0"));
    int memory_type = stoi(get_with_default(flags, "memory_type", "0")); // all disk
    if (op_type == "bulk") {
        test_blipp_bulk(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                        count, has_size);
    } else if (op_type == "lookup") {
        test_blipp_lookup(memory_type, const_cast<char *>(index_name.c_str()),
                          const_cast<char *>(key_file_path.c_str()), count, has_size, search_count, case_id, step,
                          need_sleep);
    } else if (op_type == "scan") {
        test_blipp_scan(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                        count, has_size, search_count, case_id, r_size, step);
    } else if (op_type == "insert") {
        test_insert(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                    count, has_size);
    } else if (op_type == "mix_workload") {
        hybrid_test(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                    count, has_size, case_id);
    } else if (op_type == "bulk_search_range") {
        test_bulk_search(memory_type, const_cast<char *>(index_name.c_str()), const_cast<char *>(key_file_path.c_str()),
                         count, has_size, search_count, case_id, r_size);
    }
    return 0;
}