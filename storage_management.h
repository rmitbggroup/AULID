#include<map>
#include<cstring>
#include<stack>
#include "lipp.h"
#include<vector>
#include <forward_list>

const long BLOCK_SIZE = 8192 / 2;

//typedef uint8_t bitmap_t;
//#define BITMAP_WIDTH (sizeof(bitmap_t) * 8)
//#define BITMAP_SIZE(num_items) (((num_items) + BITMAP_WIDTH - 1) / BITMAP_WIDTH)
//#define BITMAP_GET(bitmap, pos) (((bitmap)[(pos) / BITMAP_WIDTH] >> ((pos) % BITMAP_WIDTH)) & 1)
//#define BITMAP_SET(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] |= 1 << ((pos) % BITMAP_WIDTH))
//#define BITMAP_CLEAR(bitmap, pos) ((bitmap)[(pos) / BITMAP_WIDTH] &= ~bitmap_t(1 << ((pos) % BITMAP_WIDTH)))
//#define BITMAP_NEXT_1(bitmap_item) __builtin_ctz((bitmap_item))


#define v2 1

class MemoryManagement {

private:
    // 8,16,32,64
    // 0,1,2,3
    int block_info[5] = {-1, -1, -1, -1, -1};
    int used_info[5] = {0, 0, 0, 0, 0};
    int max_used[5] = {BLOCK_SIZE / 64, BLOCK_SIZE / 128, BLOCK_SIZE / 256, BLOCK_SIZE / 512, BLOCK_SIZE / 1024};
    int offset_uint[5] = {64, 128, 256, 512, 1024};
    int max_items[5] = {4, 8, 16, 32, 64};
    int total_created[5] = {0, 0, 0, 0, 0};

    class _addr {
    public:
        int block;
        int offset;
    public:
        _addr(int x, int y) : block(x), offset(y) {}
    };

    std::forward_list<_addr> p4;

    std::forward_list<_addr> p8;
//        int p8_size = 0;
    std::forward_list<_addr> p16;
//        int p16_size = 0;
    std::forward_list<_addr> p32;
//        int p32_size = 0;
    std::forward_list<_addr> p64;
//        int p64_size = 0;
    std::forward_list<int> pb;
//        int pb_size = 0;
    int current_invalid_size[6] = {0, 0, 0, 0, 0, 0};
    int MAX_INVALID = 2000;


public:
    int TOTAL_BLOCK = 0;

    void print() {
        std::cout << total_created[0] << ";" << total_created[1] << ";" << total_created[2] << ";"
                  << total_created[3] << ";" << total_created[4] << std::endl;
    }

    int get_next_addr(int type_i, int *block, int *offset) {
        total_created[type_i] += 1;
        //if (get_invalid(type_i, block, offset)) return 1;
        if (block_info[type_i] == -1 || used_info[type_i] == max_used[type_i])
            return 0;
        *block = block_info[type_i];
        *offset = offset_uint[type_i] * used_info[type_i];
        used_info[type_i] += 1;
        return 1;
    }

    void set_next_block(int type_i, int block) {
        TOTAL_BLOCK += 1;
        block_info[type_i] = block;
        used_info[type_i] = 1;
        return;
    }

    bool is_full(int type_i, int cur_num) {
        // the first item is used for item count
        return (cur_num + 1) >= max_items[type_i];
    }

    bool get_invalid(int type_i, int *block, int *offset) {
//            return false;
        if (current_invalid_size[type_i] <= 0) return false;
        current_invalid_size[type_i] -= 1;
        switch (type_i) {
            case 0: {
                auto _x = p4.front();
                *block = _x.block;
                *offset = _x.offset;
                p4.pop_front();
                break;
            }
            case 1: {
                auto _x = p8.front();
                *block = _x.block;
                *offset = _x.offset;
                p8.pop_front();
                break;
            }
            case 2: {
                auto _x = p16.front();
                *block = _x.block;
                *offset = _x.offset;
                p16.pop_front();
                break;
            }
            case 3: {
                auto _x = p32.front();
                *block = _x.block;
                *offset = _x.offset;
                p32.pop_front();
                break;
            }
            case 4: {
                auto _x = p64.front();
                *block = _x.block;
                *offset = _x.offset;
                p64.pop_front();
                break;
            }
            case 5: {
                auto _x = pb.front();
                *block = _x;
//                    *block = _x.block;
//                    *offset = _x.offset;
                pb.pop_front();
                break;
            }
            default: {
                std::cout << "WRONG... type" << std::endl;
            }
        }
        return true;
    }

    void add_invalid(int type_i, int block, int offset) {
        if (current_invalid_size[type_i] >= MAX_INVALID) return;
        current_invalid_size[type_i] += 1;
        switch (type_i) {
            case 0: {
                p4.emplace_front(block, offset);
                break;
            }
            case 1: {
                p8.emplace_front(block, offset);
                break;
            }
            case 2: {
                p16.emplace_front(block, offset);
                break;
            }
            case 3: {
                p32.emplace_front(block, offset);
                break;
            }
            case 4: {
                p64.emplace_front(block, offset);
                break;
            }
            case 5: {
                pb.emplace_front(block);
            }
            default: {
                std::cout << "WRONG... type" << std::endl;
            }
        }
    }

};

typedef struct {
    long offset; // also used as id
    int begin_items_c;
    int total_items_c;
    int data_items_c;
    int level;
} LIPPNodeMeta;

class MetaManagement {
    LIPPNodeMeta *meta = new LIPPNodeMeta[200];
    std::map<long, int> addresses;

public:

};

typedef struct {
    int is_two;
    int size;
    int build_size;//must
    int fixed;
    int number_inserts;
    int num_insert_to_data;
    int num_items; //must
    // int empty_size2; // for items
    double slope;
    long double intercept;
} NodeHeaderD;
int NodeHeaderDSize = sizeof(NodeHeaderD);

template<class T, class P>
class LIPPBTree {

    const double BUILD_LR_REMAIN = 0;

    inline int compute_gap_count(int size, bool is_firs_time = true) {
        if (size >= 1000000) return 1;
        if (size >= 100000) return 5;
        if (is_firs_time) return 10;
        else return 5;
//        return 1;
    }

    FILE *fp;

    // B-Tree related structures
    typedef struct {
        short item_count;
        int next_block_id;
        int last_block_id;
    } BLeafNodeHeader;
    const int BLeafNodeHeaderSize = sizeof(BLeafNodeHeader);

    typedef struct {
        T key;
        P value;
    } BItem;
    const int BItemSize = sizeof(BItem);
    const int MaxItemSizeInBLeafNode = (BLOCK_SIZE - BLeafNodeHeaderSize) / BItemSize;


    typedef struct {
        char tag; // see below
        union {
            struct { // 3
                T key;
                int value;
            } data;
            struct { // 2
                int block;
                int offset;
                NodeHeaderD child_head;
            } addr;
            struct { // 4-8,5-16,6-32,7-64
                int block;
                int offset;
                // todo: move to the header of actual data
                // int cur_num;
            } packed_array;
            struct { // 8
                int block1;
                int block2;
                int block3;
                int block4;
                // todo: move to the header of actual data
//                uint16_t cur_num1;
//                uint16_t cur_num2;
//                uint16_t cur_num3;
//                uint16_t cur_num4;
                T key1;
                T key2;
                T key3;
            } b_tree;
        } comp;
    } ItemD;
    int ItemDSize = sizeof(ItemD);
    int MaxItemCount = BLOCK_SIZE / ItemDSize;

    // TODO: currently, suppose the root node must be lipp node (case 2).
    //  Later, it can be different e.g., insert from scratch.
    typedef struct {
        int root_block;
        int root_offset;
        int NEXT_BLOCK;
        int NEXT_OFFSET;
        T lower_key;
        T upper_key;
        int block_id;
        NodeHeaderD root_head;
    } MetaNode;
    int MetaNodeSize = sizeof(MetaNode);

    typedef struct {
        T key;
        int value;
    } KVPair;

    MemoryManagement _mm_;

    int NEXT_BLOCK;
    int NEXT_OFFSET;
    MetaNode mb;
    LIPP<T, int> lipp_inner;

#define ALL_DISK 0
#define LEAF_DISK 1
    int MemoryType = ALL_DISK;

private:
    void write_data(void *data, long offset, int len) {
        fseek(fp, offset, SEEK_SET);
        fwrite(data, len, 1, fp);
        return;
    }

    void read_block(void *data, int block_id) {
        fseek(fp, block_id * BLOCK_SIZE, SEEK_SET);
        fread(data, BLOCK_SIZE, 1, fp);
        return;
    }

    void write_block(void *data, int block_id) {
        fseek(fp, block_id * BLOCK_SIZE, SEEK_SET);
        fwrite(data, BLOCK_SIZE, 1, fp);
        return;
    }

    void read_data(void *data, long offset, int len) {
        fseek(fp, offset, SEEK_SET);
        fread(data, len, 1, fp);
        return;
    }


public:
    //StorageManager (
    void print_mm() {
        _mm_.print();
    }

    size_t report_file_size() {
        fseek(fp, 0, SEEK_END);
        return ftell(fp);
    }

    void print_mb() {
        std::cout << mb.NEXT_BLOCK << std::endl;
    }

    int TOTAL_LEAF = 0;
    int TOTAL_BTREE = 0;

    void sys_metablock(bool update_root, NodeHeaderD head, int block = 0, int offset = 0) {
        char empty_block[BLOCK_SIZE];
        mb.NEXT_BLOCK = NEXT_BLOCK;
        mb.NEXT_OFFSET = NEXT_OFFSET;
        if (update_root) {
            mb.root_block = block;
            mb.root_offset = offset;
            memcpy(&mb.root_head, &head, NodeHeaderDSize);
        }
        memcpy(empty_block, &mb, MetaNodeSize);
        write_data(empty_block, 0, BLOCK_SIZE);
    }

    void load_metablock() {
        char data[BLOCK_SIZE];
        read_block(data, 0);
        MetaNode *_mb = (MetaNode *) (data);
        mb.NEXT_BLOCK = _mb->NEXT_BLOCK;
        mb.NEXT_OFFSET = _mb->NEXT_OFFSET;
        mb.root_block = _mb->root_block;
        mb.root_offset = _mb->root_offset;
        mb.block_id = _mb->block_id;
        mb.lower_key = _mb->lower_key;
        mb.upper_key = _mb->upper_key;
        memcpy(&mb.root_head, &(_mb->root_head), NodeHeaderDSize);
    }


    void init(char *fn, bool is_first, int memory_type) {
        if (is_first) {
            fp = fopen(fn, "wb");
            NEXT_BLOCK = 1;
            NEXT_OFFSET = 0;
            NodeHeaderD _head;
            // index.sys_metablock(false, _head);
            sys_metablock(true, _head, 1, 0);
            fclose(fp);
        }
        fp = fopen(fn, "r+b");
        load_metablock();
        MemoryType = memory_type;
    }

    void _write_node_internal(NodeHeaderD nhd, ItemD *its, int *block, int *offset, NodeHeaderD *r_head) {
        // if (BLOCK_SIZE - NEXT_OFFSET < ItemDSize) {
        //     char empty[BLOCK_SIZE - NEXT_OFFSET];
        //     write_data(empty, NEXT_BLOCK * BLOCK_SIZE + NEXT_OFFSET, BLOCK_SIZE - NEXT_OFFSET);
        //     NEXT_BLOCK += 1;
        //     NEXT_OFFSET = 0;
        // }
        // *block = NEXT_BLOCK;
        // *offset = NEXT_OFFSET;
        // // Write Header
        // long start_offset = NEXT_BLOCK * BLOCK_SIZE + NEXT_OFFSET;
        // long _start_header = start_offset;
        // write_data(&nhd, start_offset, NodeHeaderDSize);
        // start_offset += NodeHeaderDSize;

        // nhd.empty_size2 = (BLOCK_SIZE - (start_offset % BLOCK_SIZE)) % ItemDSize;
        // if (start_offset % BLOCK_SIZE == 0) {
        //     nhd.empty_size2 = 0;
        // }
        // write_data(&nhd, _start_header, NodeHeaderDSize);
        // memcpy(r_head, &nhd, NodeHeaderDSize);
        // char empty[nhd.empty_size2];
        // write_data(empty, start_offset, nhd.empty_size2);
        // start_offset += nhd.empty_size2;
        memcpy(r_head, &nhd, NodeHeaderDSize);
        {
            if (NEXT_OFFSET > 0) {
                NEXT_BLOCK += 1;
                NEXT_OFFSET = 0;
            }
        }
        long start_offset = NEXT_BLOCK * BLOCK_SIZE + NEXT_OFFSET;
        long empty_size = (BLOCK_SIZE - NEXT_OFFSET) % ItemDSize;

        // 默认想的时，如果一个完全用来存数据，则会存在开头，尾部被填充为0
        // 但是现在头部被填充为0
        if (empty_size > 0) {
            char empty[empty_size];
            write_data(empty, start_offset, empty_size);
        }
        start_offset += empty_size;
        *block = start_offset / BLOCK_SIZE;
        *offset = start_offset % BLOCK_SIZE;
        int c1 = (BLOCK_SIZE - (start_offset % BLOCK_SIZE)) / ItemDSize;
        // if (c1 < MaxItemCount) {
        if (c1 > nhd.num_items) c1 = nhd.num_items;
        write_data(its, start_offset, ItemDSize * c1);
        start_offset += ItemDSize * c1;
        // } else {
        //     c1 = 0;
        // }

        // NEXT_BLOCK = start_offset / BLOCK_SIZE;
        char data[BLOCK_SIZE];
        for (int i = c1; i < nhd.num_items;) {
            int _c = nhd.num_items - i;
            if (_c > MaxItemCount) _c = MaxItemCount;
            if (_c == MaxItemCount) {
                memcpy(data, its + i, _c * ItemDSize);
                write_data(data, start_offset, BLOCK_SIZE);
                start_offset += BLOCK_SIZE;
            } else {
                write_data(its + i, start_offset, _c * ItemDSize);
                start_offset += _c * ItemDSize;
            }
            i += _c;
        }

        NEXT_BLOCK = start_offset / BLOCK_SIZE;
        NEXT_OFFSET = start_offset % BLOCK_SIZE;
        NodeHeaderD _head;
        sys_metablock(false, _head);
    }

    typedef struct {
        int block;
        int offset;
    } Address;

    bool _search_leaf_node(char *data, int l_block, T key, T itd_key, P *r_value, int *bc, int *IPOS = nullptr,
                           int *BLOCK_ID = nullptr) {
        // char data[BLOCK_SIZE];
        read_block(data, l_block);
        *bc += 1;
        BLeafNodeHeader *blnh = (BLeafNodeHeader *) data;
        BItem *bis = (BItem *) (data + BLeafNodeHeaderSize);
        // std::cout.precision(17);
        // std::cout << key << std::endl;
        while (itd_key < key && blnh->next_block_id != -1) {
            l_block = blnh->next_block_id;
            read_block(data, l_block);
            *bc += 1;
            itd_key = bis[blnh->item_count - 1].key;
            //std::cout << "FROM TRVERSAL" << std::endl;
        }
        if (itd_key == key) {
            if (IPOS != nullptr) {
                // memcpy(DATA, data, BLOCK_SIZE);
                *IPOS = blnh->item_count - 1;
                *BLOCK_ID = l_block;
            } else *r_value = bis[blnh->item_count - 1].value;
            return true;
        }
        // binary search
        bool found = false;
        int l = 0;
        int r = blnh->item_count - 1;
        while (l <= r) {
            int mid = l + (r - l) / 2;
            if (bis[mid].key == key) {
                found = true;
                if (IPOS != nullptr) {
                    // memcpy(DATA, data, BLOCK_SIZE);
                    *IPOS = mid;
                    *BLOCK_ID = l_block;
                } else *r_value = bis[mid].value;
                return found;
            } else if (bis[mid].key < key) l = mid + 1;
            else r = mid - 1;
        }
        if (IPOS != nullptr) {
            // memcpy(DATA, data, BLOCK_SIZE);
            *BLOCK_ID = l_block;
            *IPOS = l;
        }

        return found;
    }

    bool
    _get_next_data_item_from_childv3(NodeHeaderD child_header, char *data2, bool same, int block, int offset, T key,
                                     P *r_value, int *bc, int *inner_c, int *IPOS = nullptr, int *BLOCK_ID = nullptr) {
        char data[BLOCK_SIZE];
        if (!same) {
            read_block(data, block);
            *bc += 1;
        } else memcpy(data, data2, BLOCK_SIZE);

        int last_block = block;
        int item_count = child_header.num_items;
        long _offset_ = (last_block * BLOCK_SIZE + offset) % BLOCK_SIZE;
        int c1 = (BLOCK_SIZE - _offset_) / ItemDSize;
        *inner_c += 1;


        for (int i = 0; i < item_count;) {
            if (last_block != block) {
                read_block(data, block);
                last_block = block;
                *bc += 1;
            }
            ItemD *_items = (ItemD *) (data + _offset_);
            if (c1 > (item_count - i)) {
                c1 = item_count - i;
            }
            for (int j = 0; j < c1; j++) {
                if (_items[j].tag == 3) {
                    return _search_leaf_node(data2, _items[j].comp.data.value, key,
                                             _items[j].comp.data.key, r_value, bc, IPOS, BLOCK_ID);
                } else if (_items[j].tag == 2) {
                    if (_items[j].comp.addr.block == block) memcpy(data2, data, BLOCK_SIZE);
                    return _get_next_data_item_from_childv3(_items[j].comp.addr.child_head, data2,
                                                            _items[j].comp.addr.block == block,
                                                            _items[j].comp.addr.block,
                                                            _items[j].comp.addr.offset, key, r_value,
                                                            bc, inner_c, IPOS, BLOCK_ID);
                } else if (_items[j].tag == 4 || _items[j].tag == 5 || _items[j].tag == 6 || _items[j].tag == 7 ||
                           _items[j].tag == 8) {
                    int _o_ = _items[j].comp.packed_array.offset;
                    read_block(data, _items[j].comp.packed_array.block);
                    KVPair *kvs = (KVPair *) (data + _o_);
                    *inner_c += 1;
                    *bc += 1;
                    int i_c = kvs[0].value;
                    for (int kv_i = 1; kv_i < i_c + 1; kv_i++) {
                        if (kvs[kv_i].key < key) continue;
                        return _search_leaf_node(data2, kvs[kv_i].value, key,
                                                 kvs[kv_i].key, r_value, bc, IPOS, BLOCK_ID);
                    }
                    return _search_leaf_node(data2, kvs[i_c].value, key,
                                             kvs[i_c].key, r_value, bc, IPOS, BLOCK_ID);
                } else if (_items[j].tag == 9) {
                    // read which block
                    int to_access_block = _items[j].comp.b_tree.block1;
                    if (_items[j].comp.b_tree.block2 > 0 && _items[j].comp.b_tree.key1 < key)
                        to_access_block = _items[j].comp.b_tree.block2;
                    if (_items[j].comp.b_tree.block3 > 0 && _items[j].comp.b_tree.key2 < key)
                        to_access_block = _items[j].comp.b_tree.block3;
                    if (_items[j].comp.b_tree.block4 > 0 && _items[j].comp.b_tree.key3 < key)
                        to_access_block = _items[j].comp.b_tree.block4;
                    // read block & find the leaf node address
                    read_block(data, to_access_block);
                    KVPair *kvs = (KVPair *) data;
                    *inner_c += 1;
                    *bc += 1;
                    int i_c = kvs[0].value;
                    for (int kv_i = 1; kv_i < i_c + 1; kv_i++) {
                        if (kvs[kv_i].key < key) continue;
                        return _search_leaf_node(data2, kvs[kv_i].value, key,
                                                 kvs[kv_i].key, r_value, bc, IPOS, BLOCK_ID);
                    }
                    return _search_leaf_node(data2, kvs[i_c].value, key,
                                             kvs[i_c].key, r_value, bc, IPOS, BLOCK_ID);
                }
            }
            i += c1;
            block += 1;
            _offset_ = 0;
            c1 = MaxItemCount;
        }
        return false;
    }

#define FulFill 0
#define ScanF 1

    bool _search_branch2(NodeHeaderD child_header, char *data2, bool same, bool *found_leaf, int block, int offset,
                         T key, P *r_value, int *bc, int *inner_c, int *IPOS = nullptr, int *BLOCK_ID = nullptr) {
        char data[BLOCK_SIZE];
        // obtain the next item
        int pos = predict_pos(child_header, key);
        long _offset_ = (block * BLOCK_SIZE + offset) % BLOCK_SIZE;
        int c1 = (BLOCK_SIZE - _offset_) / ItemDSize;
        ItemD *itd;
        *inner_c += 1;
        if (_offset_ == 0) {
            c1 = 0;
        }
        if (pos < c1) {
            if (!same) {
                read_block(data, block);
                *bc += 1;
//                    *inner_c += 1;
            } else memcpy(data, data2, BLOCK_SIZE);
            itd = (ItemD *) (data + _offset_ + pos * ItemDSize);
            _offset_ = pos;
        } else {

            if (c1 > 0) block += 1;

            int _b = (pos - c1) / MaxItemCount;
            block += _b;
            _offset_ = (pos - c1) % MaxItemCount;
            read_block(data, block);
            *bc += 1;

//                *inner_c += 1;

            itd = (ItemD *) (data + _offset_ * ItemDSize);

            c1 = MaxItemCount;
            if (child_header.num_items - pos < c1 - _offset_) {
                c1 = child_header.num_items - pos + _offset_;
            }
        }
        // go to different branchs based on item type
        if (itd->tag == 2) {
            if (itd->comp.addr.block == block) memcpy(data2, data, BLOCK_SIZE);
            bool found = _search_branch2(itd->comp.addr.child_head, data2,
                                         itd->comp.addr.block == block, found_leaf,
                                         itd->comp.addr.block, itd->comp.addr.offset,
                                         key, r_value, bc, inner_c, IPOS, BLOCK_ID);
            if (*found_leaf) return found;
        } else if (itd->tag == 4 || itd->tag == 5 || itd->tag == 6 || itd->tag == 7 || itd->tag == 8) {
//                if (itd->comp.packed_array.block == block) memcpy(data2, data, BLOCK_SIZE);
            bool found = search_lippb_p_v3(itd->tag, *itd, data2,
                                           false, found_leaf,
                                           itd->comp.packed_array.block,
                                           itd->comp.packed_array.offset, key, r_value, bc,
                                           inner_c, IPOS, BLOCK_ID);
            if (*found_leaf) return found;
        } else if (itd->tag == 9) {
            bool found = search_lippb_p_v3(itd->tag, *itd, data2, // read the block from itd
                                           false, found_leaf, -1, -1, key, r_value, bc,
                                           inner_c, IPOS, BLOCK_ID);
            if (*found_leaf) return found;
            return false;
        } else if (itd->tag == 3) {
            int l_block = itd->comp.data.value;
            T l_key = itd->comp.data.key;
#if ScanF
            {
                // to reduce the fetched blocks, we hope we can find the next `data' slot in the same block...
                // otherwise, we just use current found one.
                if (key > l_key) {
                    _offset_ += 1;
                    itd += 1;
                    pos += 1;
                    for (;;) {
                        if (_offset_ < c1 && pos < child_header.num_items) {
                            if (itd->tag == 3) {
                                l_block = itd->comp.data.value;
                                l_key = itd->comp.data.key;
                                break;
                            } else if (itd->tag == 4 || itd->tag == 5 || itd->tag == 6
                                       || itd->tag == 7 || itd->tag == 8 || itd->tag == 9) {
                                break;
                            }
#if FulFill
                            l_block = itd->comp.data.value;
                            l_key = itd->comp.data.key;
                            break;
#else
                            // for tag 1
                            _offset_ += 1;
                            pos += 1;
                            itd += 1;
#endif
                        } else {
                            break;
                        }
                    }
                }
            }
#endif
            *found_leaf = true;
            return _search_leaf_node(data2, l_block, key, l_key, r_value, bc, IPOS, BLOCK_ID);
        } else if (itd->tag == 1) {
#if FulFill
            *found_leaf = true;
            return _search_leaf_node(data2, itd->comp.data.value, key, itd->comp.data.key, r_value, bc, IPOS, BLOCK_ID);
#else
            _offset_ += 1;
            itd += 1;
            pos += 1;
            for (;;) {
                if (_offset_ < c1 && pos < child_header.num_items) {
                    if (itd->tag == 3) {
                        *found_leaf = true;
                        return _search_leaf_node(data2, itd->comp.data.value, key, itd->comp.data.key, r_value, bc,
                                                 IPOS, BLOCK_ID);
                    } else if (itd->tag == 2) {
                        *found_leaf = true;
                        if (itd->comp.addr.block == block) memcpy(data2, data, BLOCK_SIZE);
                        return _get_next_data_item_from_childv3(itd->comp.addr.child_head, data2,
                                                                itd->comp.addr.block == block, itd->comp.addr.block,
                                                                itd->comp.addr.offset, key, r_value, bc, inner_c, IPOS,
                                                                BLOCK_ID);
                    } else if (itd->tag == 4 || itd->tag == 5 || itd->tag == 6 || itd->tag == 7 || itd->tag == 8) {
                        *found_leaf = true;
                        *bc += 1;
                        *inner_c += 1;
                        read_block(data2, itd->comp.packed_array.block);
                        KVPair *kvs = (KVPair *) (data2 + itd->comp.packed_array.offset);
                        int i_c = kvs[0].value;
                        for (int kv_i = 1; kv_i < i_c + 1; kv_i++) {
                            if (kvs[kv_i].key < key) continue;
                            return _search_leaf_node(data2, kvs[kv_i].value, key,
                                                     kvs[kv_i].key, r_value, bc, IPOS, BLOCK_ID);
                        }
                        return _search_leaf_node(data2, kvs[i_c].value, key,
                                                 kvs[i_c].key, r_value, bc, IPOS, BLOCK_ID);
                    } else if (itd->tag == 9) {
                        // read which block
                        *found_leaf = true;
                        *bc += 1;
                        *inner_c += 1;

                        int to_access_block = itd->comp.b_tree.block1;
                        if (itd->comp.b_tree.block2 > 0 && itd->comp.b_tree.key1 < key)
                            to_access_block = itd->comp.b_tree.block2;
                        if (itd->comp.b_tree.block3 > 0 && itd->comp.b_tree.key2 < key)
                            to_access_block = itd->comp.b_tree.block3;
                        if (itd->comp.b_tree.block4 > 0 && itd->comp.b_tree.key3 < key)
                            to_access_block = itd->comp.b_tree.block4;
                        // read block & find the leaf node address
                        read_block(data2, to_access_block);
                        KVPair *kvs = (KVPair *) data2;
                        int i_c = kvs[0].value;
                        for (int kv_i = 1; kv_i < i_c + 1; kv_i++) {
                            if (kvs[kv_i].key < key) continue;
                            return _search_leaf_node(data2, kvs[kv_i].value, key,
                                                     kvs[kv_i].key, r_value, bc, IPOS, BLOCK_ID);
                        }
                        return _search_leaf_node(data2, kvs[i_c].value, key,
                                                 kvs[i_c].key, r_value, bc, IPOS, BLOCK_ID);
                        //std::cout << "WRONG AREA AT _search_branch2 - tag-" << itd->tag << std::endl;
                    }
                    _offset_ += 1;
                    pos += 1;
                    itd += 1;
                } else if (pos < child_header.num_items) { // next block
                    //std::cout << "READ FROM EMTPY..." << std::endl;
                    block += 1;
                    read_block(data, block);
                    itd = (ItemD *) data;
                    *bc += 1;
                    _offset_ = 0;
                    c1 = MaxItemCount;
                    if (child_header.num_items - pos < MaxItemCount) {
                        c1 = child_header.num_items - pos;
                    }
                    continue;
                } else { //go up
                    *found_leaf = false;
                    return false;
                }
            }
#endif
        }
    }

    bool search_lippb_p_v3(int type_i, ItemD _item_, char *data2, bool same, bool *found_leaf, int block, int offset,
                           T key, P *r_value, int *bc, int *inner_c, int *IPOS = nullptr, int *BLOCK_ID = nullptr) {
#if v2
        if (key >= mb.lower_key) {
            // printf("last...\n");
            *found_leaf = true;
            return _search_leaf_node(data2, mb.block_id, key, mb.upper_key, r_value, bc, IPOS, BLOCK_ID);
        }
#endif
        switch (type_i) {
            case 1: {
                *found_leaf = false;
                std::cout << "WRONG CALL AT search_lippb_p_v3  case 1";
                return false;
                break;
            }

            case 2: {
                return _search_branch2(_item_.comp.addr.child_head, data2,
                                       same, found_leaf, block, offset, key,
                                       r_value, bc, inner_c, IPOS, BLOCK_ID);
            }

            case 3: {
                std::cout << "WRONG CALL AT search_lippb_p_v3 case 3";
                return false;
                break; // data
            }
            case 4: //8
            case 5: //16
            case 6: //32
            case 7: //64
            case 8: {
//                    std::cout << "FROM PACKED" << std::endl;
                *found_leaf = true;
                read_block(data2, _item_.comp.packed_array.block);
                *bc += 1;
                *inner_c += 1;
                KVPair *kvs = (KVPair *) (data2 + _item_.comp.packed_array.offset);
                int i_c = kvs[0].value;
                for (int kv_i = 1; kv_i < i_c + 1; kv_i++) {
                    if (kvs[kv_i].key < key) continue;
                    return _search_leaf_node(data2, kvs[kv_i].value, key,
                                             kvs[kv_i].key, r_value, bc, IPOS, BLOCK_ID);
                }
                return _search_leaf_node(data2, kvs[i_c].value, key,
                                         kvs[i_c].key, r_value, bc, IPOS, BLOCK_ID);
            }
            case 9: {
//                    std::cout << "FROM B+-tree NODE" << std::endl;
                *found_leaf = true;
                int to_access_block = _item_.comp.b_tree.block1;
                if (_item_.comp.b_tree.block2 > 0 && _item_.comp.b_tree.key1 < key)
                    to_access_block = _item_.comp.b_tree.block2;
                if (_item_.comp.b_tree.block3 > 0 && _item_.comp.b_tree.key2 < key)
                    to_access_block = _item_.comp.b_tree.block3;
                if (_item_.comp.b_tree.block4 > 0 && _item_.comp.b_tree.key3 < key)
                    to_access_block = _item_.comp.b_tree.block4;
                // read block & find the leaf node address
                read_block(data2, to_access_block);
                *bc += 1;
                *inner_c += 1;

                KVPair *kvs = (KVPair *) data2;
                int i_c = kvs[0].value;
                for (int kv_i = 1; kv_i < i_c + 1; kv_i++) {
                    if (kvs[kv_i].key < key) continue;
                    return _search_leaf_node(data2, kvs[kv_i].value, key,
                                             kvs[kv_i].key, r_value, bc, IPOS, BLOCK_ID);
                }
                return _search_leaf_node(data2, kvs[i_c].value, key,
                                         kvs[i_c].key, r_value, bc, IPOS, BLOCK_ID);
            }
            default:
                std::cout << "WRONG BRANCH AT search_lippb_p_v3 default" << std::endl;
                return false;
        }

    }

    bool lipp_search_all_disk(T key, P *r_value, int *bc, int *inner_c) {
        bool found_leaf = false;
        char data[BLOCK_SIZE];
        ItemD item;
        item.comp.addr.child_head = mb.root_head;
        item.comp.addr.block = mb.root_block;
        item.comp.addr.offset = mb.root_offset;
        bool found = search_lippb_p_v3(2, item, data, false, &found_leaf, mb.root_block, mb.root_offset, key, r_value,
                                       bc, inner_c);
        if (!found)
            return _search_leaf_node(data, mb.block_id, key, mb.upper_key, r_value, bc, nullptr, nullptr);
        if (!found_leaf) {
            printf("bug...\n");
        }
        return found;
    }

    bool lipp_search_leaf_disk(char *data, T key, P *r_value, int *bc, int *inner_c, int *IPOS = nullptr,
                               int *BLOCK_ID = nullptr) {
        bool found_leaf = false;
//            char data[BLOCK_SIZE];
#if v2
        if (key >= mb.lower_key) {
            // printf("last...\n");
            found_leaf = true;
//                char data[BLOCK_SIZE];
            return _search_leaf_node(data, mb.block_id, key, mb.upper_key, r_value, bc, IPOS, BLOCK_ID);;
        }
        std::vector<std::pair<T, int>> results;
        int found = lipp_inner.range_query_len(results, key, 1);
        found_leaf = found == 1;
        if (!found_leaf)
            return _search_leaf_node(data, mb.block_id, key, mb.upper_key, r_value, bc, IPOS, BLOCK_ID);
        return _search_leaf_node(data, results[0].second, key, results[0].first, r_value, bc, IPOS, BLOCK_ID);
#endif
    }

    bool lippb_search_entry(T key, P *r_value, int *bc, int *inner_c) {
        // search_lippb(bool *found_leaf, int block, int offset, T key, P *r_value, int *bc) {
        if (MemoryType == ALL_DISK) {
            return lipp_search_all_disk(key, r_value, bc, inner_c);
        } else {
            char data[BLOCK_SIZE];
            return lipp_search_leaf_disk(data, key, r_value, bc, inner_c);
        }
    }

    void lippb_scan_entry(T key, T *keys, int *bc, int len, int *inner_c) {
        char start_block[BLOCK_SIZE];
        bool found_leaf = false;
        int start_pos = 0;
        int block_id = 0;
        ItemD item;
        item.comp.addr.child_head = mb.root_head;
        item.comp.addr.block = mb.root_block;
        item.comp.addr.offset = mb.root_offset;
        if (MemoryType == ALL_DISK)
            search_lippb_p_v3(2, item, start_block, false, &found_leaf, mb.root_block, mb.root_offset, key, nullptr, bc,
                              inner_c, &start_pos, &block_id);
        else
            lipp_search_leaf_disk(start_block, key, nullptr, bc, inner_c, &start_pos, &block_id);
        BLeafNodeHeader *blfnh = (BLeafNodeHeader *) start_block;
        // int count = blfnh->item_count;
        BItem *bis = (BItem *) (start_block + BLeafNodeHeaderSize);
        // std::cout << i << ","<<  start_pos <<"," <<  blfnh->item_count << std::endl;
        // std::cout << *bc << "," << start_pos << "," << blfnh->item_count << std::endl;
        // std::cout << len << std::endl;
        int i = 0;
        for (; i < len;) {
            for (; start_pos < blfnh->item_count && i < len; start_pos++, i++) {
                keys[i] = bis[start_pos].key;
            }
            if (i < len && blfnh->next_block_id != -1) {
                read_block(start_block, blfnh->next_block_id);
                *bc += 1;
                start_pos = 0;
            }
        }
    }

    int predict_pos(NodeHeaderD nhd, T key) {
        double v = nhd.slope * static_cast<long double>(key) + nhd.intercept;
        int pos = 0;
        if (v > std::numeric_limits<int>::max() / 2) {
            pos = nhd.num_items - 1;
        } else if (v < 0) {
            pos = 0;
        } else {
            pos = std::min(nhd.num_items - 1, static_cast<int>(v));
        }
        return pos;
    }


    void print_block() {
        std::cout << TOTAL_BTREE << "," << TOTAL_LEAF << "," << _mm_.TOTAL_BLOCK << std::endl;
    }

    void bulk_load_entry(T *_keys, P *_values, int _size, int *r_block = nullptr, int *r_offset = nullptr) {

        TOTAL_BTREE = TOTAL_LEAF = 0;
        _mm_.TOTAL_BLOCK = 0;
        int INSERT_SIZE_LEAF = int(MaxItemSizeInBLeafNode * 0.8);
        int leaf_count = _size / INSERT_SIZE_LEAF;
        int last_item_count = _size % INSERT_SIZE_LEAF;
        if (last_item_count > 0) leaf_count += 1;
#if v2
        T *last_keys = new T[leaf_count - 1]; // scan forward
        int *addresses = new int[leaf_count - 1];
        std::vector<std::pair<T, int>> data_lipp;
#else
        T last_keys[leaf_count]; // scan forward
        int addresses[leaf_count];
#endif
        char data[BLOCK_SIZE];
        if (NEXT_OFFSET > 0) {
            char empty[BLOCK_SIZE - NEXT_OFFSET];
            write_data(empty, NEXT_BLOCK * BLOCK_SIZE + NEXT_OFFSET, BLOCK_SIZE - NEXT_OFFSET);
            NEXT_BLOCK += 1;
            NEXT_OFFSET = 0;
        }
//            char *data = new char[BLOCK_SIZE];
        for (int i = 0; i < leaf_count; i++) {
            BLeafNodeHeader blnh;
            blnh.item_count = (last_item_count > 0 && i == leaf_count - 1 ? last_item_count : INSERT_SIZE_LEAF);
            blnh.next_block_id = (i == leaf_count - 1 ? -1 : NEXT_BLOCK + 1);
            blnh.last_block_id = (i == 0 ? -1 : NEXT_BLOCK - 1);
            memcpy(data, &blnh, BLeafNodeHeaderSize);
            BItem _bis[blnh.item_count];
            for (int j = 0; j < blnh.item_count; j++) {
                _bis[j].key = _keys[i * INSERT_SIZE_LEAF + j];
                _bis[j].value = _values[i * INSERT_SIZE_LEAF + j];
            }
            memcpy(data + BLeafNodeHeaderSize, _bis, BItemSize * blnh.item_count);
            write_block(data, NEXT_BLOCK);
#if v2
            if (i < leaf_count - 1) {
                last_keys[i] = _bis[blnh.item_count - 1].key;
                addresses[i] = NEXT_BLOCK;
                data_lipp.push_back({_bis[blnh.item_count - 1].key, NEXT_BLOCK});
            } else {
                mb.block_id = NEXT_BLOCK;
                mb.upper_key = _bis[blnh.item_count - 1].key;
                mb.lower_key = _bis[0].key;
            }
#else
            last_keys[i] = _bis[blnh.item_count - 1].key;
            addresses[i] = NEXT_BLOCK;
#endif
            NEXT_BLOCK += 1;
        }
        NodeHeaderD _head;
        sys_metablock(false, _head);
        printf("leaf_count:%d\n", leaf_count);
#if v2
        int LC = leaf_count - 1;
#else
        int LC = leaf_count;
#endif
        if (MemoryType == ALL_DISK) {
//                bulk_load_disk_lipp(last_keys, addresses, LC);
            long total_level = 0;
            int _b, _o;
            NodeHeaderD nhd;
            bulk_load_disk_lippv2(last_keys, addresses, 0, LC, 1, &total_level, true, &_b, &_o, &nhd);
        } else
            lipp_inner.bulk_load(data_lipp.data(), data_lipp.size());
        delete[]last_keys;
        delete[]addresses;
    }


    int WRITE_LENs[6] = {64, 128, 256, 512, 1024, BLOCK_SIZE};

    void write_packed_array(int type_i, int cur_num, int block, int offset, KVPair *kvs, bool is_virtal_one = false) {
        char data[BLOCK_SIZE];
        KVPair kv;
        kv.value = cur_num;
        memcpy(data, &kv, sizeof(KVPair));
        memcpy(data + sizeof(KVPair), kvs, cur_num * sizeof(KVPair));
        long t_offset = block * BLOCK_SIZE + offset;
        write_data(data, t_offset, WRITE_LENs[type_i]);
        return;
    }

    int bulk_load_disk_lippv2(T *_keys, int *_values, int begin, int end, int level, long *TOTAL_LEVEL,
                              bool is_first_bulk = true, int *r_block = nullptr, int *r_offset = nullptr,
                              NodeHeaderD *r_head = nullptr) {
        T *keys = _keys + begin;
        int *values = _values + begin;
        const int size = end - begin;
        const int BUILD_GAP_CNT = compute_gap_count(size, is_first_bulk);
        NodeHeaderD node;
        node.is_two = 0;
        node.build_size = size;
        node.size = size;
        node.fixed = 0;
        node.number_inserts = node.num_insert_to_data = 0;
        int DATA_L3 = 0;
//            long TOTAL_LEVEL = 0;
        if (level == 1) { // build start from the root node
            *TOTAL_LEVEL = 0;
        }
        // call FMCD
        {
            const int L = size * static_cast<int>(BUILD_GAP_CNT + 1);
            int i = 0;
            int D = 1;
            double Ut = (static_cast<long double>(keys[size - 1 - D]) - static_cast<long double>(keys[D])) /
                        (static_cast<double>(L - 2)) + 1e-6;
            while (i < size - 1 - D) {
                while (i + D < size && keys[i + D] - keys[i] >= Ut) {
                    i++;
                }
                if (i + D >= size) {
                    break;
                }
                D = D + 1;
                if (D * 3 > size) break;
                Ut = (static_cast<long double>(keys[size - 1 - D]) - static_cast<long double>(keys[D])) /
                     (static_cast<double>(L - 2)) + 1e-6;
            }
            std::cout << "CF:" << D << std::endl;
            if (D * 3 <= size) {

                node.slope = 1.0 / Ut;
                node.intercept = (L - node.slope * (static_cast<long double>(keys[size - 1 - D]) +
                                                    static_cast<long double>(keys[D]))) / 2;

                node.num_items = L;
            } else {
                int mid1_pos = (size - 1) / 3;
                int mid2_pos = (size - 1) * 2 / 3;

                const long double mid1_key = (static_cast<long double>(keys[mid1_pos]) +
                                              static_cast<long double>(keys[mid1_pos + 1])) / 2;
                const long double mid2_key = (static_cast<long double>(keys[mid2_pos]) +
                                              static_cast<long double>(keys[mid2_pos + 1])) / 2;

                node.num_items = size * static_cast<int>(BUILD_GAP_CNT + 1);
                const double mid1_target =
                        mid1_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;
                const double mid2_target =
                        mid2_pos * static_cast<int>(BUILD_GAP_CNT + 1) + static_cast<int>(BUILD_GAP_CNT + 1) / 2;

                node.slope = (mid2_target - mid1_target) / (mid2_key - mid1_key);
                node.intercept = mid1_target - node.slope * mid1_key;
            }
        }
        const int lr_remains = static_cast<int>(size * BUILD_LR_REMAIN);
        node.intercept += lr_remains;
        node.num_items += lr_remains * 2;

        if (size > 1e6) {
            node.fixed = 1;
        }
        ItemD *items = new ItemD[node.num_items];
        for (int i = 0; i < node.num_items; i++) {
            items[i].tag = 1;
        }

        for (int item_i = predict_pos(node, keys[0]), offset = 0; offset < size;) {

            int next = offset + 1, next_i = -1;
            while (next < size) {
                next_i = predict_pos(node, keys[next]);
                if (next_i == item_i) {
                    next++;
                } else {
                    break;
                }
            }

#if FulFill
            {
                        for (int f_i = item_i - 1; f_i >=0; f_i--) {
                            if (items[f_i].tag == 1) {
                                items[f_i].comp.data.key = keys[offset];
                                items[f_i].comp.data.value = values[offset];
                            } else {
                                break;
                            }
                        }
                    }
#endif
            int insertion_count = next - offset;
            if (insertion_count == 1) {
                if (level >= 3) DATA_L3 += 1;
//#if CalAvg
                *TOTAL_LEVEL += level * 1;
//#endif
                items[item_i].tag = 3;
                items[item_i].comp.data.key = keys[offset];
                items[item_i].comp.data.value = values[offset];

            } else if (insertion_count < 64) {
//#if CalAvg
                *TOTAL_LEVEL += insertion_count * (level + 1);
//#endif
                if (level >= 2) DATA_L3 += insertion_count;
                int l_block;
                int l_offset;

                KVPair *kvs = new KVPair[insertion_count];
                for (int c_i = 0; c_i < insertion_count; c_i++) {
                    kvs[c_i].key = keys[offset + c_i];
                    kvs[c_i].value = values[offset + c_i];
                }


                int type_i;
                if (insertion_count < 4) {
                    type_i = 0;
                } else if (insertion_count < 8) {
                    type_i = 1;
                } else if (insertion_count < 16) {
                    type_i = 2;
                } else if (insertion_count < 32) {
                    type_i = 3;
                } else {
                    type_i = 4;
                }
                items[item_i].tag = type_i + 4;
                if (_mm_.get_next_addr(type_i, &l_block, &l_offset) == 0) {
                    if (NEXT_OFFSET > 0) {
                        NEXT_BLOCK += 1;
                        NEXT_OFFSET = 0;
                    }
                    _mm_.set_next_block(type_i, NEXT_BLOCK);
                    l_block = NEXT_BLOCK;
                    l_offset = NEXT_OFFSET;
                    // after that, we keep this block for i8
                    NEXT_BLOCK += 1;
                }
                write_packed_array(type_i, insertion_count, l_block, l_offset, kvs);
                items[item_i].comp.packed_array.block = l_block;
                items[item_i].comp.packed_array.offset = l_offset;
                delete[] kvs;
            } else if (insertion_count <
                       (BLOCK_SIZE / sizeof(KVPair) - 1) * 4) { // 1, because the first one is item count
//#if CalAvg
                *TOTAL_LEVEL += insertion_count * (level + 1);
//#endif
                if (level >= 3) DATA_L3 += insertion_count;
                std::cout << "BUILD B+-tree NODe" << insertion_count << std::endl;
                // determine how many blocks we need
                int ITEM_IN_A_BLOCK = int(BLOCK_SIZE / sizeof(KVPair)) - 1;
                int needed_block = insertion_count / ITEM_IN_A_BLOCK;
                int _r_ = insertion_count % ITEM_IN_A_BLOCK;
                if (_r_ > 0) needed_block += 1;
                items[item_i].tag = 9;
                items[item_i].comp.b_tree.block1 = -1;
                items[item_i].comp.b_tree.block2 = -1;
                items[item_i].comp.b_tree.block3 = -1;
                items[item_i].comp.b_tree.block4 = -1;
                if (NEXT_OFFSET > 0) {
                    NEXT_BLOCK += 1;
                    NEXT_OFFSET = 0;
                }
                TOTAL_BTREE += 1;
                // write each block
                int _start_ = offset;
                for (int b_i = 0; b_i < needed_block; b_i++) {
                    int _end_ = (b_i == needed_block - 1) ? offset + insertion_count : _start_ + ITEM_IN_A_BLOCK;
                    int ac = _end_ - _start_;
                    KVPair *kvs = new KVPair[ac];
                    for (int c_i = _start_; c_i < _end_; c_i++) {
                        kvs[c_i - _start_].key = keys[c_i];
                        kvs[c_i - _start_].value = values[c_i];
                    }
                    write_packed_array(5, ac, NEXT_BLOCK, 0, kvs);
                    switch (b_i) {
                        case 0: {
                            items[item_i].comp.b_tree.block1 = NEXT_BLOCK;
                            items[item_i].comp.b_tree.key1 = kvs[ac - 1].key;
                            break;
                        }
                        case 1: {
                            items[item_i].comp.b_tree.block2 = NEXT_BLOCK;
                            items[item_i].comp.b_tree.key2 = kvs[ac - 1].key;
                            break;
                        }
                        case 2: {
                            items[item_i].comp.b_tree.block3 = NEXT_BLOCK;
                            items[item_i].comp.b_tree.key3 = kvs[ac - 1].key;
                            break;
                        }
                        case 3: {
                            items[item_i].comp.b_tree.block4 = NEXT_BLOCK;
                            break;
                        }
                        default:
                            std::cout << "ERROR BLOCK BRANCH at BULK" << std::endl;
                    }
                    NEXT_BLOCK += 1;
                    _start_ += ac;
                    delete[] kvs;
                }
            } else {
                std::cout << "LARGER THAN 4 BLOCKS" << std::endl;
                items[item_i].tag = 2;
                NodeHeaderD c_nhd;
                // call algo again and insert
                int c_block;
                int c_offset;
                int c_l3 = bulk_load_disk_lippv2(_keys, _values, begin + offset, begin + next, level + 1,
                                                 TOTAL_LEVEL, is_first_bulk, &c_block, &c_offset, &c_nhd);
                items[item_i].comp.addr.block = c_block;
                items[item_i].comp.addr.offset = c_offset;
                memcpy(&(items[item_i].comp.addr.child_head), &c_nhd, NodeHeaderDSize);
                DATA_L3 += c_l3;
            }

            if (next >= size) {
                break;
            } else {
                item_i = next_i;
                offset = next;
            }
        }
        node.num_insert_to_data = DATA_L3;
        int cl_block;
        int cl_offset;
        NodeHeaderD _head;
        _write_node_internal(node, items, &cl_block, &cl_offset, &_head);
        if (level == 1) {// update the root node in main memory
            std::cout << "avg:" << *TOTAL_LEVEL / double(end) << std::endl;
            *r_block = cl_block;
            *r_offset = cl_offset;
            sys_metablock(true, _head, cl_block, cl_offset);
            memcpy(r_head, &_head, NodeHeaderDSize);
        } else {
            *r_block = cl_block;
            *r_offset = cl_offset;
            memcpy(r_head, &_head, NodeHeaderDSize);
        }
        std::cout << "file size" << report_file_size() << std::endl;
        delete[]items;
        return DATA_L3;
    }

#define CalAvg 1

    void scan_and_destory_treev3(NodeHeaderD _root, int block, int offset, T *keys, int *values) {
        typedef std::pair<int, NodeHeaderD> Segment;
        typedef std::pair<int, int> NodeAddress;
        std::stack<Segment> s;
        std::stack<NodeAddress> addr;

        s.push(Segment(0, _root));
        addr.push(NodeAddress(block, offset));
        char *data = new char[BLOCK_SIZE];
        char *data2 = new char[BLOCK_SIZE];
        while (!s.empty()) {
            int begin = s.top().first;
            NodeHeaderD node = s.top().second;
            const int SHOULD_END_POS = begin + node.size;
            s.pop();
            int _block = addr.top().first;
            int _offset = addr.top().second;
            addr.pop();
            long _offset_ = (_block * BLOCK_SIZE + _offset) % BLOCK_SIZE;
            int c1 = (BLOCK_SIZE - _offset_) / ItemDSize;

            // if (_offset_ == 0) {
            //     _block += 1;
            // }

            for (int i = 0; i < node.num_items;) {
                read_block(data, _block);
                ItemD *_items = (ItemD *) (data + _offset_);
                if (c1 > (node.num_items - i)) {
                    c1 = node.num_items - i;
                }
                for (int j = 0; j < c1; j++) {
                    if (_items[j].tag == 3) {
                        keys[begin] = _items[j].comp.data.key;
                        values[begin] = _items[j].comp.data.value;
                        begin++;
                    } else if (_items[j].tag == 2) {
                        int _b = _items[j].comp.addr.block;
                        int _o = _items[j].comp.addr.offset;
                        NodeHeaderD _n;
                        memcpy(&_n, &_items[j].comp.addr.child_head, NodeHeaderDSize);
                        s.push(Segment(begin, _n));
                        addr.push(NodeAddress(_b, _o));
                        begin += _n.size;
                    } else if (_items[j].tag == 4 || _items[j].tag == 5 || _items[j].tag == 6 || _items[j].tag == 7) {
                        int _b = _items[j].comp.packed_array.block;
                        int _o = _items[j].comp.packed_array.offset;
                        read_block(data2, _b);
                        KVPair *kvs = (KVPair *) (data2 + _o);
                        int _ic = kvs[0].value;
                        for (int _i_ = 1; _i_ < _ic + 1; _i_++) {
                            keys[begin] = kvs[_i_].key;
                            values[begin] = kvs[_i_].value;
                            begin++;
                        }
                    } else if (_items[j].tag == 8) {
                        if (_items[j].comp.b_tree.block1 > 0) {
                            read_block(data2, _items[j].comp.b_tree.block1);
                            KVPair *kvs = (KVPair *) (data2);
                            int _ic = kvs[0].value;
                            for (int _i_ = 1; _i_ < _ic + 1; _i_++) {
                                keys[begin] = kvs[_i_].key;
                                values[begin] = kvs[_i_].value;
                                begin++;
                            }
                        }
                        if (_items[j].comp.b_tree.block2 > 0) {
                            read_block(data2, _items[j].comp.b_tree.block2);
                            KVPair *kvs = (KVPair *) (data2);
                            int _ic = kvs[0].value;
                            for (int _i_ = 1; _i_ < _ic + 1; _i_++) {
                                keys[begin] = kvs[_i_].key;
                                values[begin] = kvs[_i_].value;
                                begin++;
                            }
                        }
                        if (_items[j].comp.b_tree.block3 > 0) {
                            read_block(data2, _items[j].comp.b_tree.block3);
                            KVPair *kvs = (KVPair *) (data2);
                            int _ic = kvs[0].value;
                            for (int _i_ = 1; _i_ < _ic + 1; _i_++) {
                                keys[begin] = kvs[_i_].key;
                                values[begin] = kvs[_i_].value;
                                begin++;
                            }
                        }
                        if (_items[j].comp.b_tree.block4 > 0) {
                            read_block(data2, _items[j].comp.b_tree.block4);
                            KVPair *kvs = (KVPair *) (data2);
                            int _ic = kvs[0].value;
                            for (int _i_ = 1; _i_ < _ic + 1; _i_++) {
                                keys[begin] = kvs[_i_].key;
                                values[begin] = kvs[_i_].value;
                                begin++;
                            }
                        }

                    }
                }

                i += c1;
                _block += 1;
                _offset_ = 0;
                c1 = MaxItemCount;
            }
        }
        delete[] data;
        delete[] data2;
    }


    void count_root() {
        int num_items = 0;
        int num_nodes = 0;
        int insert_count = 0;
        count_one_node(mb.root_block, mb.root_offset, &num_items, &num_nodes, &insert_count);
        std::cout << num_items << std::endl;
        std::cout << num_nodes << std::endl;
        std::cout << insert_count << std::endl;
        // std::cout << count_one_node(mb.root_block, mb.root_offset) << std::endl;
    }

    void count_one_node(int block, int offset, int *NumItems, int *NumNodes, int *InsertCount) {
        char data[BLOCK_SIZE];
        char data2[BLOCK_SIZE];
        read_block(data2, block);
        //read one node header
        int _block = block;

        // int count  = 0;
        // int node_count = 1;
        NodeHeaderD *nhd = (NodeHeaderD *) (data2 + offset);
        long _offset_ = (_block * BLOCK_SIZE + offset) % BLOCK_SIZE;
        int c1 = (BLOCK_SIZE - _offset_) / ItemDSize;

        // if (_offset_ == 0) {
        //     _block += 1;
        // }

        // std::cout<<nhd->num_items<<std::endl;
        *NumItems += nhd->num_items;
        *NumNodes += 1;
        for (int i = 0; i < nhd->num_items;) {
            read_block(data, _block);
            ItemD *_items = (ItemD *) (data + _offset_);
            if (c1 > (nhd->num_items - i)) {
                c1 = nhd->num_items - i;
            }
            for (int j = 0; j < c1; j++) {
                if (_items[j].tag == 2) {
                    count_one_node(_items[j].comp.addr.block, _items[j].comp.addr.offset, NumItems, NumNodes,
                                   InsertCount);
                } else if (_items[j].tag == 3) {
                    *InsertCount += 1;
                }
            }
            i += c1;
            _block += 1;
            _offset_ = 0;
            c1 = MaxItemCount;
        }
    }

    bool insert_lippb_entry(const T &key, const P &value, int *bc, int *in_c) {
        // find the leaf node to insert new key
        char start_block[BLOCK_SIZE];
        bool found_leaf = false;
        int start_pos = 0;
        int block_id = 0;
        *in_c = 0;
        if (MemoryType == ALL_DISK) {
            ItemD item;
            item.comp.addr.child_head = mb.root_head;
            item.comp.addr.block = mb.root_block;
            item.comp.addr.offset = mb.root_offset;
            search_lippb_p_v3(2, item, start_block, false, &found_leaf, mb.root_block, mb.root_offset, key,
                              nullptr, bc, in_c, &start_pos, &block_id);
            if (!found_leaf) {
                _search_leaf_node(start_block, mb.block_id, key, mb.upper_key, nullptr, bc, &start_pos,
                                  &block_id); // should not happen
            }
        } else {
            lipp_search_leaf_disk(start_block, key, nullptr, bc, in_c, &start_pos, &block_id);
        }
        BLeafNodeHeader *blfnh = (BLeafNodeHeader *) start_block;
        BItem *bis = (BItem *) (start_block + BLeafNodeHeaderSize);
        while (start_pos < blfnh->item_count) {
            if (bis[start_pos].key < key) start_pos++;
            else break;
        }
        if (blfnh->item_count == MaxItemSizeInBLeafNode) {
            if (NEXT_OFFSET > 0) {
                char empty[BLOCK_SIZE - NEXT_OFFSET];
                write_data(empty, NEXT_BLOCK * BLOCK_SIZE + NEXT_OFFSET, BLOCK_SIZE - NEXT_OFFSET);
                NEXT_BLOCK += 1;
                NEXT_OFFSET = 0;
            }
            TOTAL_LEAF += 1;
            BItem temp[MaxItemSizeInBLeafNode + 1];
            memcpy(temp, bis, MaxItemSizeInBLeafNode * BItemSize);
            for (int i = blfnh->item_count - 1; i >= start_pos; i--) {
                temp[i + 1] = temp[i];
            }
            temp[start_pos].key = key;
            temp[start_pos].value = value;

            char data[BLOCK_SIZE];
            BLeafNodeHeader *_blnh = (BLeafNodeHeader *) data;

            _blnh->item_count = (MaxItemSizeInBLeafNode + 1) / 2;
            _blnh->next_block_id = block_id;
            _blnh->last_block_id = blfnh->last_block_id;
            blfnh->item_count = MaxItemSizeInBLeafNode + 1 - _blnh->item_count;
            blfnh->last_block_id = NEXT_BLOCK;

            BItem *_bis = (BItem *) (data + BLeafNodeHeaderSize);
            // to avoid `update' operation here, we keep the larger part in the original node.
            for (int i = 0; i < _blnh->item_count; i++) {
                _bis[i] = temp[i];
            }

            for (int i = _blnh->item_count, j = 0; i < MaxItemSizeInBLeafNode + 1; i++, j++) {
                bis[j] = temp[i];
            }
#if v2
            if (block_id == mb.block_id) { // the last segment => update the lower and upper key
                mb.lower_key = bis[0].key;
                mb.upper_key = bis[blfnh->item_count - 1].key;
            }
#endif
            write_block(start_block, block_id);
            write_block(data, NEXT_BLOCK);


            int last_block = _blnh->last_block_id;
            T new_key = _bis[_blnh->item_count - 1].key;
            int new_addr = NEXT_BLOCK;
            if (last_block != -1) {
                read_block(data, last_block);
                BLeafNodeHeader *_last_blnh = (BLeafNodeHeader *) data;
                _last_blnh->next_block_id = NEXT_BLOCK;
                write_block(data, last_block);
            }
            NEXT_BLOCK += 1;
            // sys_metablock(false);
            if (MemoryType == ALL_DISK)
                return insertv3(new_key, new_addr);
            else {
                lipp_inner.insert(new_key, new_addr);
                return true;
            }
        } else {
            for (int i = blfnh->item_count - 1; i >= start_pos; i--) {
                bis[i + 1] = bis[i];
            }
            bis[start_pos].key = key;
            bis[start_pos].value = value;
            blfnh->item_count += 1;
            // check
            if (block_id == mb.block_id && (start_pos == 0 || start_pos == blfnh->item_count - 1)) {
                mb.lower_key = bis[0].key;
                mb.upper_key = bis[blfnh->item_count - 1].key;
                // sys_metablock(false);
            }
            write_block(start_block, block_id);
            return true;
        }

    }

    void get_next_item_for_node2(NodeHeaderD node_header, int block, int offset, char *parent_block_data,
                                 int *parent_block, int *item_offset, T key) {
        int pos = predict_pos(node_header, key);
        long _offset_ = (block * BLOCK_SIZE + offset) % BLOCK_SIZE;
        int c1 = (BLOCK_SIZE - _offset_) / ItemDSize;
        if (_offset_ == 0)
            c1 = 0;
        if (pos < c1) {
            read_block(parent_block_data, block);
            *parent_block = block;
            *item_offset = _offset_ + pos * ItemDSize;
        } else {
            int _block_ = c1 > 0 ? block + 1 : block;
            _block_ += (pos - c1) / MaxItemCount;
            read_block(parent_block_data, _block_);
            *parent_block = _block_;
            *item_offset = ((pos - c1) % MaxItemCount) * ItemDSize;
        }
        return;
    }

    // caller makes sure having enough space...
    // the first item in kvs is the number of item in kvs
    void merge_an_item(KVPair *new_kvs, KVPair *kvs, int cur_num, T key, int value) {
        bool is_added = false;
        for (int n_i = 0, o_i = 1; n_i < cur_num + 1 && o_i < cur_num + 1; n_i++) {
            if (!is_added && kvs[o_i].key > key) {
                new_kvs[n_i].key = key;
                new_kvs[n_i].value = value;
                is_added = true;
                continue;
            }
            new_kvs[n_i].key = kvs[o_i].key;
            new_kvs[n_i].value = kvs[o_i].value;
            o_i += 1;
        }
        if (!is_added) {
            new_kvs[cur_num].key = key;
            new_kvs[cur_num].value = value;
        }
        return;
    }

    // currently, we suppose the root node is lipp node type
    bool insertv3(const T &key, const int &value) {
//            if (value == 21539) {
//                int x = 0;
//            }
//            char *parent_block_data = new char[BLOCK_SIZE];
//            int parent_block = -1;
//            int item_offset = -1;
        char *data = new char[BLOCK_SIZE];
        ItemD *header_item[5]; // if more than 5 lipp nodes we pass, not meaningful
        int header_blocks[5];
        int header_offsets[5];
        char *PointerToBlock[5];
        int INSERT_INTO_L3 = 0;
        int header_i = 0; //header_i is root
        int valid_header_ic = 0;
        bool last_lipp_is_gen = false;

        {
//                memcpy(header_traj + header_i, &mb.root_head, NodeHeaderDSize);
            ItemD *itd_h = new ItemD;
            itd_h->comp.addr.child_head = mb.root_head;
            itd_h->tag = 2;
            itd_h->comp.addr.block = mb.root_block;
            itd_h->comp.addr.offset = mb.root_offset;
//                header_traj[header_i] = &mb.root_head;
            header_item[header_i] = itd_h;
            header_blocks[header_i] = mb.root_block;
            header_offsets[header_i] = mb.root_offset;
            header_i += 1;
            valid_header_ic += 1;
        }

        // fetch item
        ItemD *item_d;
        {
            int parent_block = -1;
            int item_offset = -1;
            char *parent_block_data = new char[BLOCK_SIZE];
            get_next_item_for_node2(mb.root_head, mb.root_block,
                                    mb.root_offset, parent_block_data,
                                    &parent_block, &item_offset, key);
            PointerToBlock[header_i] = parent_block_data;
            header_blocks[header_i] = parent_block;
            header_offsets[header_i] = item_offset;
        }
        item_d = (ItemD *) (PointerToBlock[header_i] + header_offsets[header_i]);
        for (;;) {
            bool is_insert_suc = false;
            switch (item_d->tag) {
                case 1: // meet empty & do insertion
                {
                    item_d->tag = 3;
                    item_d->comp.data.key = key;
                    item_d->comp.data.value = value;
                    write_block(PointerToBlock[header_i], header_blocks[header_i]);
                    is_insert_suc = true;
                    if (header_i >= 3) {
                        INSERT_INTO_L3 += 1;
                    }

                    break;
//                        delete[] parent_block_data;
//                        delete[] data;
//                        return true;
                }
                case 2: // meet lipp node again
                {
                    {
//                            header_traj[header_i] = &(item_d->comp.addr.child_head);
                        header_item[header_i] = item_d;
                        valid_header_ic += 1;

                        header_i += 1;
                        int parent_block = -1;
                        int item_offset = -1;
                        char *parent_block_data = new char[BLOCK_SIZE];
                        get_next_item_for_node2(item_d->comp.addr.child_head, item_d->comp.addr.block,
                                                item_d->comp.addr.offset, parent_block_data,
                                                &parent_block, &item_offset, key);
                        PointerToBlock[header_i] = parent_block_data;
                        header_blocks[header_i] = parent_block;
                        header_offsets[header_i] = item_offset;
                    }
                    item_d = (ItemD *) (PointerToBlock[header_i] + header_offsets[header_i]);
                    break;
                }
                case 3: // data & meet conflict
                {
                    KVPair *kvs = new KVPair[2];
                    kvs[0].key = item_d->comp.data.key;
                    kvs[0].value = item_d->comp.data.value;
                    kvs[1].key = key;
                    kvs[1].value = value;
                    if (item_d->comp.data.key > key) {
                        kvs[1].key = item_d->comp.data.key;
                        kvs[1].value = item_d->comp.data.value;
                        kvs[0].key = key;
                        kvs[0].value = value;
                    }
                    int l_block, l_offset;
                    if (_mm_.get_next_addr(0, &l_block, &l_offset) == 0) {
                        if (NEXT_OFFSET > 0) {
                            NEXT_BLOCK += 1;
                            NEXT_OFFSET = 0;
                        }
                        _mm_.set_next_block(0, NEXT_BLOCK);
                        l_block = NEXT_BLOCK;
                        l_offset = NEXT_OFFSET;
                        NEXT_BLOCK += 1;
                    }
                    write_packed_array(0, 2, l_block, l_offset, kvs);
                    // update the item, and we also need to update its parent?
                    item_d->tag = 4;
                    item_d->comp.packed_array.block = l_block;
                    item_d->comp.packed_array.offset = l_offset;
                    write_block(PointerToBlock[header_i], header_blocks[header_i]);
                    is_insert_suc = true;
                    if (header_i >= 2) {
                        INSERT_INTO_L3 += 2;
                    }
                    delete[] kvs;
                    break;
//                        delete[] parent_block_data;
//                        delete[] data;
//                        return true;
                }
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:// packed array & do insertion & upgrade if needed
                {
                    if (header_i >= 2) {
                        INSERT_INTO_L3 += 1;
                    }

                    int _block_ = item_d->comp.packed_array.block;
//                        if (_block_ == 41751) {
//                            std::cout << 0 << std::endl;
//                        }
                    int _offset_ = item_d->comp.packed_array.offset;
                    int type_i = item_d->tag - 4;
                    read_block(data, _block_);
                    KVPair *kvs = (KVPair *) (data + _offset_);
                    int cur_num = kvs[0].value;

                    // merge into the new array
                    KVPair *new_kvs = new KVPair[cur_num + 1];
                    merge_an_item(new_kvs, kvs, cur_num, key, value);

                    // directly insert into the orginal node
                    if (!_mm_.is_full(type_i, cur_num)) {
                        KVPair _kv_;
                        _kv_.value = cur_num + 1;
                        memcpy(data + _offset_, &_kv_, sizeof(KVPair));
                        memcpy(data + _offset_ + sizeof(KVPair),
                               new_kvs, sizeof(KVPair) * (cur_num + 1));
                        write_block(data, _block_);
                        delete[] new_kvs;
                        is_insert_suc = true;
                        break;
//                            delete[] parent_block_data;
//                            delete[] data;
//                            return true;
                    }

                    // re-write the current node into a new node
                    //_mm_.add_invalid(type_i, _block_, _offset_);
                    if (item_d->tag < 8) { // upgrade to a higher packed array
                        int l_block, l_offset;
                        if (_mm_.get_next_addr(type_i + 1, &l_block, &l_offset) == 0) {
                            if (NEXT_OFFSET > 0) {
                                NEXT_BLOCK += 1;
                                NEXT_OFFSET = 0;
                            }
                            _mm_.set_next_block(type_i + 1, NEXT_BLOCK);
                            l_block = NEXT_BLOCK;
                            l_offset = NEXT_OFFSET;
                            NEXT_BLOCK += 1;
                        }
                        write_packed_array(type_i + 1, cur_num + 1, l_block,
                                           l_offset, new_kvs);
                        item_d->comp.packed_array.block = l_block;
                        item_d->comp.packed_array.offset = l_offset;
                    } else { // upgrade to b+-tree node
                        if (NEXT_OFFSET > 0) {
                            NEXT_BLOCK += 1;
                            NEXT_OFFSET = 0;
                        }
                        TOTAL_BTREE += 1;
                        item_d->comp.b_tree.block1 = NEXT_BLOCK;
                        item_d->comp.b_tree.block2 = -1;
                        item_d->comp.b_tree.block3 = -1;
                        item_d->comp.b_tree.block4 = -1;
                        item_d->comp.b_tree.key1 = new_kvs[cur_num].key;
                        write_packed_array(5, cur_num + 1, NEXT_BLOCK, 0, new_kvs);
                        NEXT_BLOCK += 1;
                    }
                    delete[] new_kvs;
                    // update the parent node
                    item_d->tag = item_d->tag + 1;
                    write_block(PointerToBlock[header_i], header_blocks[header_i]);
//                        delete[] parent_block_data;
//                        delete[] data;
                    is_insert_suc = true;
                    break;
//                        return true;
                }
                case 9: // b_tree node & do insertion & convert to lipp node if needed
                {
                    int MAX_ITEM = BLOCK_SIZE / sizeof(KVPair) - 1;
                    int to_insert_block = item_d->comp.b_tree.block1;
                    int to_insert_index = 1;
                    if (item_d->comp.b_tree.block2 > 0 && item_d->comp.b_tree.key1 < key) {
                        to_insert_block = item_d->comp.b_tree.block2;
                        to_insert_index = 2;
                    }
                    if (item_d->comp.b_tree.block3 > 0 && item_d->comp.b_tree.key2 < key) {
                        to_insert_block = item_d->comp.b_tree.block3;
                        to_insert_index = 3;
                    }
                    if (item_d->comp.b_tree.block4 > 0 && item_d->comp.b_tree.key3 < key) {
                        to_insert_block = item_d->comp.b_tree.block4;
                        to_insert_index = 4;
                    }
                    int used_block = 1 + (item_d->comp.b_tree.block2 > 0)
                                     + (item_d->comp.b_tree.block3 > 0)
                                     + (item_d->comp.b_tree.block4 > 0);
                    read_block(data, to_insert_block);
                    KVPair *kvs = (KVPair *) data;
                    int i_c = kvs[0].value;

                    if (i_c < MAX_ITEM) { // insert into that block
                        int o_i = i_c;
                        for (; o_i > 0; o_i--) {
                            if (kvs[o_i].key < key) break;
                            kvs[o_i + 1] = kvs[o_i];
                        }
                        kvs[o_i + 1].key = key;
                        kvs[o_i + 1].value = value;
                        kvs[0].value = i_c + 1;
                        write_block(data, to_insert_block);
//                            delete[] parent_block_data;
//                            delete[] data;
//                            return true;
                        is_insert_suc = true;
                        if (header_i >= 2) {
                            INSERT_INTO_L3 += 1;
                        }
                        break;
                    } else if (i_c == MAX_ITEM && used_block < 4) { // do split
                        if (header_i >= 2) {
                            INSERT_INTO_L3 += 1;
                        }
                        // merge
                        KVPair *new_kvs = new KVPair[i_c + 1];
                        merge_an_item(new_kvs, kvs, i_c, key, value);
                        int i_c_l = (i_c + 1) / 2;
                        int i_c_r = (i_c + 1) - i_c_l;
                        write_packed_array(5, i_c_l, to_insert_block, 0, new_kvs);
                        if (NEXT_OFFSET > 0) {
                            NEXT_BLOCK += 1;
                            NEXT_OFFSET = 0;
                        }
                        TOTAL_BTREE += 1;
                        write_packed_array(5, i_c_r, NEXT_BLOCK, 0, new_kvs + i_c_l);

                        if (to_insert_index == 4) {
                            std::cout << "INSERT ERROR ... AT insertv3 - case 8";
//                                delete[] parent_block_data;
//                                delete[] data;
//                                return false;
                            is_insert_suc = true;
                            break;
                        } else if (to_insert_index == 3) {
                            item_d->comp.b_tree.block4 = NEXT_BLOCK;
                            item_d->comp.b_tree.key3 = new_kvs[i_c_l - 1].key;
                        } else if (to_insert_index == 2) {
                            if (item_d->comp.b_tree.block3 > 0) {
                                item_d->comp.b_tree.block4 = item_d->comp.b_tree.block3;
                                item_d->comp.b_tree.key3 = item_d->comp.b_tree.key2;
                            }
                            item_d->comp.b_tree.block3 = NEXT_BLOCK;
                            item_d->comp.b_tree.key2 = new_kvs[i_c_l - 1].key;
                        } else {
                            if (item_d->comp.b_tree.block3 > 0) {
                                item_d->comp.b_tree.block4 = item_d->comp.b_tree.block3;
                                item_d->comp.b_tree.key3 = item_d->comp.b_tree.key2;
                            }
                            if (item_d->comp.b_tree.block2 > 0) {
                                item_d->comp.b_tree.block3 = item_d->comp.b_tree.block2;
                                item_d->comp.b_tree.key2 = item_d->comp.b_tree.key1;
                            }
                            item_d->comp.b_tree.block2 = NEXT_BLOCK;
                            item_d->comp.b_tree.key1 = new_kvs[i_c_l - 1].key;
                        }
                        write_block(PointerToBlock[header_i], header_blocks[header_i]);
                        NEXT_BLOCK += 1;
                        delete[] new_kvs;
//                            delete[] parent_block_data;
//                            delete[] data;
//                            return true;
                        is_insert_suc = true;
                        break;
                    } else { // fetch all items in 4 blocks and build lipp node
                        //bulk_load_disk_lipp()
                        std::cout << "SPECIAL CASE UPGRADE B+-TREE NODE" << std::endl;
                        T *_keys_ = new T[MAX_ITEM * 4];
                        int *_values_ = new int[MAX_ITEM * 4];
                        read_block(data, item_d->comp.b_tree.block1);
                        KVPair *kvs = (KVPair *) data;
                        int i_c = kvs[0].value;
                        int g_i = 0;
                        bool is_added = false;
                        for (int l_i = 1; l_i < i_c + 1; g_i += 1) {
                            if (!is_added && kvs[l_i].key > key) {
                                _keys_[g_i] = key;
                                _values_[g_i] = value;
                                is_added = true;
                                continue;
                            }
                            _keys_[g_i] = kvs[l_i].key;
                            _values_[g_i] = kvs[l_i].value;
                            l_i++;
                        }

                        read_block(data, item_d->comp.b_tree.block2);
                        kvs = (KVPair *) data;
                        i_c = kvs[0].value;
                        for (int l_i = 1; l_i < i_c + 1; g_i += 1) {
                            if (!is_added && kvs[l_i].key > key) {
                                _keys_[g_i] = key;
                                _values_[g_i] = value;
                                is_added = true;
                                continue;
                            }
                            _keys_[g_i] = kvs[l_i].key;
                            _values_[g_i] = kvs[l_i].value;
                            l_i++;
                        }

                        read_block(data, item_d->comp.b_tree.block3);
                        kvs = (KVPair *) data;
                        i_c = kvs[0].value;
                        for (int l_i = 1; l_i < i_c + 1; g_i += 1) {
                            if (!is_added && kvs[l_i].key > key) {
                                _keys_[g_i] = key;
                                _values_[g_i] = value;
                                is_added = true;
                                continue;
                            }
                            _keys_[g_i] = kvs[l_i].key;
                            _values_[g_i] = kvs[l_i].value;
                            l_i++;
                        }

                        read_block(data, item_d->comp.b_tree.block4);
                        kvs = (KVPair *) data;
                        i_c = kvs[0].value;
                        for (int l_i = 1; l_i < i_c + 1; g_i += 1) {
                            if (!is_added && kvs[l_i].key > key) {
                                _keys_[g_i] = key;
                                _values_[g_i] = value;
                                is_added = true;
                                continue;
                            }
                            _keys_[g_i] = kvs[l_i].key;
                            _values_[g_i] = kvs[l_i].value;
                            l_i++;
                        }


                        int _b_ = 0;
                        int _o_ = 0;
                        NodeHeaderD _header_;
//                            int DATA_SLOT_COUNT = bulk_load_disk_lipp(_keys_, _values_, g_i,
//                                                false,&_b_,&_o_, &_header_);
                        long total_levels = 0;
                        int data_l3 = bulk_load_disk_lippv2(_keys_, _values_, 0, g_i, header_i + 1, &total_levels,
                                                            false, &_b_, &_o_, &_header_);
                        item_d->tag = 2;
                        item_d->comp.addr.block = _b_;
                        item_d->comp.addr.offset = _o_;
                        memcpy(&(item_d->comp.addr.child_head), &_header_, NodeHeaderDSize);
//                            header_traj[header_i] = &(item_d->comp.addr.child_head);
                        header_item[header_i] = item_d;
//                            header_blocks[header_i] = _b_;
//                            header_offsets[he]
                        valid_header_ic += 1;
                        last_lipp_is_gen = true;
                        // leave update the header type 2 to post proccessing
//                            write_block(PointerToBlock[header_i], header_blocks[header_i]);
//                            if (header_i == 1) {
                        INSERT_INTO_L3 += data_l3;
//                            } else if (header_i > 1){
//                                INSERT_INTO_L3 += 1;
//                            }
                        delete[] _keys_;
                        delete[] _values_;
//                            delete[] parent_block_data;
//                            delete[] data;
//                            return true;
                        is_insert_suc = true;
                        break;
                    }
                }
                default: {
//                        delete[] parent_block_data;
//                        delete[] data;
                    std::cout << "WRONG BRANCH ... AT insertv3" << std::endl;
//                        return false;
                    is_insert_suc = true;
                    break;
                }
            }
            if (is_insert_suc) {
                break;
            }
        }

        //Update Stats
//            for (int i = 1; i < valid_header_ic; i++) {
//                header_traj[i]->num_insert_to_data += INSERT_INTO_L3;
//                header_traj[i]->number_inserts += 1;
//                header_traj[i]->size += 1;
//                if (i > 0) { //update header, we won't update it in the above for-loop
//                    write_block(PointerToBlock[i], header_blocks[header_i]);
//                }
//            }

        // post processing
        mb.root_head.num_insert_to_data += INSERT_INTO_L3;
        mb.root_head.number_inserts += 1;
        mb.root_head.size += 1;
        for (int i = 1; i < header_i + 1; i++) {
            if (i < valid_header_ic - 1 || (!last_lipp_is_gen && i == valid_header_ic - 1)) {
                (header_item[i]->comp.addr.child_head).num_insert_to_data += INSERT_INTO_L3;
                (header_item[i]->comp.addr.child_head).number_inserts += 1;
                (header_item[i]->comp.addr.child_head).size += 1;
//                    write_block(PointerToBlock[i], header_blocks[i]);
            } else if (last_lipp_is_gen && i == valid_header_ic - 1) {
                // todo: check is = or +=
                (header_item[i]->comp.addr.child_head).num_insert_to_data += INSERT_INTO_L3; // we need to record how many in l3
            }
//                delete[] PointerToBlock[i];
        }

        //SMO check
        int update_idx = valid_header_ic - 1;
#define SMO 1
#if SMO
        for (int i = 0; i < valid_header_ic; i++) {
            NodeHeaderD nhd = (i == 0) ? mb.root_head : (header_item[i]->comp.addr.child_head);
            int start_block_ = (i == 0) ? mb.root_block : (header_item[i]->comp.addr).block;
            int start_offset_ = (i == 0) ? mb.root_offset : (header_item[i]->comp.addr).offset;
            bool need_rebuild = 20 * nhd.num_insert_to_data >= nhd.size && nhd.size >= 1.5 * nhd.build_size;
            if (need_rebuild) {
                const int ESIZE = nhd.size;
                int old_l3 = nhd.num_insert_to_data;
                T *keys = new T[ESIZE];
                int *values = new int[ESIZE];
                scan_and_destory_treev3(nhd, start_block_, start_offset_, keys, values);
                int _b = 0;
                int _o = 0;
                NodeHeaderD _header;
//                    int data_count = bulk_load_disk_lipp(keys, values, ESIZE, false, &_b, &_o, &_header);
                long total_levels = 0;
                int d_l3 = bulk_load_disk_lippv2(keys, values, 0, ESIZE,
                                                 (i + 1), &total_levels, false,
                                                 &_b, &_o, &_header);
                // printf("re:%d;%d\n",i,_nhd.size);
                delete[] keys;
                delete[] values;
                // update the stats at header
                update_idx = i;
                if (i == 0) {
                    std::cout << "update root" << std::endl;
                    _header.num_insert_to_data = d_l3;
                    sys_metablock(true, _header, _b, _o);
                } else if (i < 3) {
                    std::cout << "UPDATE Node with:" << ESIZE << std::endl;
                    // write it at last, out of loop
                    for (int inner_i = 1; inner_i < i; inner_i++) {
                        header_item[i]->comp.addr.child_head.num_insert_to_data += (d_l3 - old_l3);
                    }
                    mb.root_head.num_insert_to_data += (d_l3 - old_l3);
                    _header.num_insert_to_data = d_l3;
                    memcpy(&(header_item[i]->comp.addr.child_head), &_header, NodeHeaderDSize);
                    header_item[i]->comp.addr.block = _b;
                    header_item[i]->comp.addr.offset = _o;
//                        if (i == 1) {
//                            mb.root_head.num_insert_to_data += (data_count - old_l3);
//                        }
//                        else { // do not neeed to update the upper level, all are larger than or equal to 3
//                            std::cout << "level - :" << i << std::endl;
//                        }
                } else {
                    std::cout << "UPDATE Node with l>3" << ESIZE << std::endl;
                    _header.num_insert_to_data = d_l3;
                    memcpy(&(header_item[i]->comp.addr.child_head), &_header, NodeHeaderDSize);
                    header_item[i]->comp.addr.block = _b;
                    header_item[i]->comp.addr.offset = _o;
                }
                break;
            }

        }
#endif
        for (int i = 1; i < update_idx + 1; i++) {
            write_block(PointerToBlock[i], header_blocks[i]);
        }
        for (int i = 1; i < header_i + 1; i++) {
            delete[]PointerToBlock[i];
        }
        delete[] data;
        return true;
    }
};