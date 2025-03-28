///////////////////////////////////////////////////
//
// ██████╗  █████╗ ██████╗ ██╗██████╗     ███████╗██╗   ██╗ ██████╗        
// ██╔══██╗██╔══██╗██╔══██╗██║██╔══██╗    ██╔════╝██║   ██║██╔═══██╗       
// ██████╔╝███████║██████╔╝██║██║  ██║    ███████╗██║   ██║██║   ██║       
// ██╔══██╗██╔══██║██╔═══╝ ██║██║  ██║    ╚════██║╚██╗ ██╔╝██║   ██║       
// ██║  ██║██║  ██║██║     ██║██████╔╝    ███████║ ╚████╔╝ ╚██████╔╝       
// ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝     ╚══════╝  ╚═══╝   ╚═════╝    
//  
//  MIT License
//  
//  Copyright (c) 2025 Severi Suominen
//  
//  Permission is hereby granted, free of charge, to use, copy, modify, merge,
//  publish, distribute, sublicense, and/or sell copies of this software, 
//  provided that the above copyright notice and this permission notice appear 
//  in all copies or substantial portions of the software.
//  
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
//  FITNESS FOR A PARTICULAR PURPOSE, AND NONINFRINGEMENT.
//
//  GitHub: https://github.com/SeveriSuominen
//
////////////////////////////////////////////////////

#pragma once



#include <cstdint>
#include <array>
#include <queue>
#include <cassert>

#include "libmorton/morton3D.h"
#include "glm/glm.hpp"

namespace rapid_svo
{
    /////////////////////////////
    // (1) 
    // to reduce complexity we want to only work in unsigned even integer space
    //
    // (2) 
    // just row major ordering node index
    //
    // (3) 
    // not a fan of using this kind of multi-line macros but its absolutely necessary
    // for code readibility, tolerance for any extra overhead is zero and we also dont
    // want to rely on compiler inlining as its never 100% guaranteed 

    #define SVO_NODE_RELATION_MATH_IMPL()\
        auto node_transformed = node_position;\
        node_transformed <<= 1;\
        const auto dist_vector = voxel_transformed - node_transformed;\
        const auto cell_vector = dist_vector / node_extent;\
        auto index = (cell_vector[0] << 2) + (cell_vector[1] << 1) + cell_vector[2];\
        auto next_node_position = node_position + (cell_vector * (component_type)(node_extent) / (component_type)2);\
        const uint8_t child_bit = (1 << index);\
        const auto exist = (node_->_mask & child_bit) != 0;
        
    ///////////////////////////////

    enum bit_width
    {
        morton_16b,
        morton_32b
    };

    namespace util 
    {
        template<typename T, size_t N>
        struct vargs_fixed
        {
        private:
            T (_args)[N];
        public:
            template<typename... ARGS, typename = std::enable_if_t<sizeof...(ARGS)==N>>
            constexpr vargs_fixed(ARGS&&... args):_args(args...) {}
            T& operator [] (const size_t i) {return _args[i];}
            constexpr size_t size() {return N;}
            T* ptr() {return &_args[0];}
        };

        template<uint8_t beg_val, uint8_t size_val, typename T>
        inline static void bset(T& n_, T value_) 
        requires (std::is_integral_v<T>) {
            assert(beg_val+size_val < sizeof(T)*8);
            const T width_ = size_val;
            const T cmask_ = ~(((1<<(width_))-1)<<beg_val);
            const T vmask_ = (value_&((1<<width_)-1))<<beg_val;
            n_ = (n_&cmask_)|vmask_;
        }

        template<uint8_t bel_val, uint8_t size_val, typename T>
        inline static void bget(T n_, T& value_) 
        requires (std::is_integral_v<T>) {
            assert(bel_val+size_val < sizeof(T)*8);
            value_ = (n_ >> bel_val) & ((1 << (size_val))-1);
        }

        template<typename T>
        inline static constexpr T log2(T val_)
        requires (std::is_integral_v<T>) {
            T res_ = 0;
            while(val_ > 1) { val_ >>= 1; ++res_; }
            return res_;
        }

        template<typename T>
        inline constexpr static T min(T a_, T b_)
        requires (std::is_integral_v<T>) {
            return b_ ^((a_^b_)&-(a_< b_));
        }

        template<typename T>
        inline constexpr static T max(T a_, T b_)
        requires (std::is_integral_v<T>) {
            return a_ ^((a_^b_)&-(a_< b_));
        }

        template<typename T>
        inline constexpr static T pad_size_to_min_alignment(T size_, T min_align_) 
        requires (std::is_integral_v<T>) {
            return (size_ + min_align_ - 1) & ~(min_align_ - 1);
        }
    }

    
    template<bit_width BIT_WIDTH>
    struct morton_util
    {
        inline static constexpr uint32_t AXIS_MAX = BIT_WIDTH == morton_16b ? 32 : 1024;

        using morton_type = std::conditional_t<BIT_WIDTH == morton_16b, uint16_t, uint32_t>;
        
        using component_type   = std::conditional_t<BIT_WIDTH == morton_16b, uint8_t, uint16_t>;
        
        using signed_component_type = std::conditional_t<BIT_WIDTH == morton_16b, int8_t, int16_t>;

        inline static void pos_to_morton(morton_type& morton, const component_type* in_) {
            morton = libmorton::m3D_e_sLUT<morton_type>(in_[0],in_[1],in_[2]);
        }

        inline static void morton_to_pos(morton_type morton, component_type* out_) {
            libmorton::m3D_d_sLUT<morton_type>(morton, out_[0],out_[1],out_[2]);
        }

        morton_util() = delete;
    };

    ///////////////////////////////
    // 16b: axis bit width 5, free to use total leftover bits 1, inc. extent 32^3
    // 32b: axis bit width 10, free to use total leftover bits 2, inc. extent inclusive 1024^3
    // 64b: axis bit width 21, free to use total leftover bits 1, inc. extent 2'097'151^3
    ///////////////////////////////
    template<typename payload_t, bit_width BIT_WIDTH,
    typename = std::enable_if_t<
        std::is_default_constructible_v <payload_t> &&
        std::is_copy_constructible_v    <payload_t> &&
        std::is_copy_assignable_v       <payload_t> &&
        std::is_move_constructible_v    <payload_t> &&
        std::is_move_assignable_v       <payload_t>>>
    struct spatial
    {
        using size_type = std::conditional_t<BIT_WIDTH == morton_16b, uint16_t, uint32_t>;

        size_type 
        _morton{};

        payload_t
        _data{};
        
        spatial() = default;

        spatial(const spatial& other_):
            _morton(other_._morton),
            _data(other_._data) 
        {};

        spatial(spatial&& other_):
            _morton(other_._morton),
            _data(std::move(other_._data)) 
        { 
            other_._data = {}; 
        };

        spatial(const morton_util<BIT_WIDTH>::component_type* in_, payload_t& data): 
            _data(data) 
        {
            encode_position(in_);
        }

        spatial(const morton_util<BIT_WIDTH>::component_type* in_, payload_t&& data): 
            _data(std::move(data)) 
        {
            encode_position(in_);
        }

        spatial& operator=(const spatial& other_) 
        {
            _morton = other_._morton;
            _data = other_._data;

            return *this;
        };

        spatial& operator=(spatial&& other_) 
        {
            if(this != &other_) {
                _morton = other_._morton;
                _data = std::move(other_._data);
                other_._data = {};
            }
            return *this;
        };
        
        spatial& 
        encode_position(const morton_util<BIT_WIDTH>::component_type* in_) {
            //assert((*in_ < morton_util<BIT_WIDTH>::AXIS_MAX).all());
            morton_util<BIT_WIDTH>::pos_to_morton(_morton, in_);
            return *this;
        }

        spatial& 
        decode_position(morton_util<BIT_WIDTH>::component_type* out_) {
            morton_util<BIT_WIDTH>::morton_to_pos(_morton, out_);
            return *this; 
        }

        [[nodiscard]] 
        constexpr uint8_t morton_size() const {
            return sizeof(size_type);
        }
    };

    template<typename T>
    using spatial_16b = spatial<T, morton_16b>;

    template<typename T>
    using spatial_32b = spatial<T, morton_32b>;
    ///////////////////////////////

    struct basic_voxel_format
    {   
        #define SVO_VFORMAT_VOXEL_IDX 0,15
        #define SVO_VFORMAT_STATE_BIT 15,1
        #define SVO_VFORMAT_TYPE_INFO 16,16
        #define SVO_VFORMAT_USER_DATA 32,32

        using pack_type = uint64_t;
        pack_type 
        _packed{};
        
        //basic_voxel_format& 
        //set_index(uint8_t* in_) {
        //    util::bset<SVO_VFORMAT_VOXEL_IDX>
        //    (_packed, static_cast<pack_type>(
        //        libmorton::m3D_e_sLUT<uint16_t>(in_[0],in_[1],in_[2])));
        //    return *this;
        //}

        //basic_voxel_format& 
        //get_index(uint8_t* out_) {
        //    pack_type val; 
        //    util::bget<SVO_VFORMAT_VOXEL_IDX>
        //    (_packed, val);
        //    libmorton::m3D_d_sLUT<uint16_t>(val, out_[0],out_[1],out_[2]);
        //    return *this; 
        //}

        basic_voxel_format& 
        set_state_bit(bool in_) {
            util::bset<SVO_VFORMAT_STATE_BIT>
            (_packed, static_cast<pack_type>(in_));
            return *this;
        }

        basic_voxel_format& 
        get_state_bit(bool& out_) {
            pack_type val; 
            util::bget<SVO_VFORMAT_VOXEL_IDX>
            (_packed, val);
            out_ = val;
            return *this;
        }

        basic_voxel_format& 
        set_type_info(uint16_t in_) {
            util::bset<SVO_VFORMAT_TYPE_INFO>
            (_packed, static_cast<pack_type>(in_));
            return *this;
        }

        basic_voxel_format& 
        get_type_info(uint16_t& out_) {
            pack_type val; 
            util::bget<SVO_VFORMAT_TYPE_INFO>
            (_packed, val);
            out_ = val;
            return *this;
        }

        basic_voxel_format& 
        set_user_data(uint32_t in_) {
            util::bset<SVO_VFORMAT_USER_DATA>
            (_packed, static_cast<pack_type>(in_));
            return *this;
        }

        basic_voxel_format& 
        get_user_data(uint32_t& out_) {
            pack_type val; 
            util::bget<SVO_VFORMAT_USER_DATA>
            (_packed, val);
            out_ = val;
            return *this;
        }
    };

    struct node_format
    {
        uint32_t
        _block_index{};
        
        uint8_t 
        _depth{};
        
        uint8_t 
        _mask{};

        [[nodiscard]] 
        bool has_children() const {
            return _mask & 0xFF;
        }
    };

    template<typename T>
    struct mem_pool 
    {
        std::vector<std::array<T, 8>>
        _blocks{};

        std::queue<uint32_t>
        _free{};

        uint32_t acquire_next_index()
        {
            if(_free.size() > 0) {
                return _free.front();
            } else {
                return _blocks.size(); 
            }
        }

        uint32_t alloc()
        {
            if(_free.size() > 0) {
                uint32_t index = _free.front();
                _free.pop();
                return index;
            } else {
                _blocks.emplace_back();
                return _blocks.size()-1; 
            }
        }

        void dealloc(uint32_t index)
        {
            _free.push(index);
        }
    };
    
    struct details_info
    {
        bool 
        _discard_overflow = false;

        std::array<uint32_t, 3> 
        _limit_max_bounds{}; 
    };

    template<bit_width BIT_WIDTH, typename FORMAT_T = basic_voxel_format, details_info DETAILS = {},
    typename = std::enable_if_t<
        std::is_default_constructible_v <FORMAT_T> &&
        std::is_copy_constructible_v    <FORMAT_T> &&
        std::is_copy_assignable_v       <FORMAT_T> &&
        std::is_move_constructible_v    <FORMAT_T> &&
        std::is_move_assignable_v       <FORMAT_T>>>
    class tree
    {   
    public:

        uint64_t byte_size()
        {
            auto size = sizeof(*this);
            size += get_node_blocks_count()  * sizeof(node_format) * 8;
            size += get_voxel_blocks_count() * sizeof(FORMAT_T) * 8;
            return size;
        }

        ////////////////////////
        // 16b: 8^5  == 32^3
        // 32b: 8^10 == 1024^3
        ////////////////////////
        inline static constexpr uint32_t 
        SPACE_ABSOLUTE_MAX_DEPTH = BIT_WIDTH == morton_16b ? 5 : 10; 

        inline static constexpr uint32_t 
        ABSOLUTE_AXIS_WIDTH = 1 << SPACE_ABSOLUTE_MAX_DEPTH;

        inline static constexpr std::array<uint16_t, 3>
        BOUNDS = {
        DETAILS._limit_max_bounds[0] > 0 ? util::min(DETAILS._limit_max_bounds[0], ABSOLUTE_AXIS_WIDTH) : ABSOLUTE_AXIS_WIDTH,
        DETAILS._limit_max_bounds[1] > 0 ? util::min(DETAILS._limit_max_bounds[1], ABSOLUTE_AXIS_WIDTH) : ABSOLUTE_AXIS_WIDTH,
        DETAILS._limit_max_bounds[2] > 0 ? util::min(DETAILS._limit_max_bounds[2], ABSOLUTE_AXIS_WIDTH) : ABSOLUTE_AXIS_WIDTH};
        
        inline static constexpr uint32_t calc_bounds_max_depth() {
            constexpr uint32_t res_ = util::log2(util::max(util::max(BOUNDS[0], BOUNDS[1]), BOUNDS[2])); 
            return util::max(res_, 3u);
        }

        inline static constexpr uint32_t 
        MAX_DEPTH = calc_bounds_max_depth(); 
        
        inline static constexpr uint32_t 
        AXIS_WIDTH = 1 << MAX_DEPTH;
       
        inline static constexpr uint32_t 
        MAX_SIZE = 
            AXIS_WIDTH*
            AXIS_WIDTH*
            AXIS_WIDTH; 
        ////////////////////////
        
        [[nodiscard]] 
        static constexpr bit_width get_type() { return BIT_WIDTH; }

        using voxel_format = FORMAT_T;

        using spatial_node = spatial<node_format, BIT_WIDTH>;

        using spatial_voxel = spatial<voxel_format, BIT_WIDTH>;
        
        using component_type = morton_util<BIT_WIDTH>::component_type;

        using signed_component_type = morton_util<BIT_WIDTH>::signed_component_type;
        
        using vector_type = glm::vec<3, component_type>;
    
    private:

        node_format
        _root_node{};

        mem_pool<node_format>
        _node_pool{};

        mem_pool<voxel_format>
        _voxel_pool{};

    public:

        tree()
        {
            _root_node = {};
            _root_node._depth = 0;
            _root_node._block_index = _node_pool.alloc();
        }

        [[nodiscard]] 
        uint32_t get_node_blocks_count() const
        {
            return _node_pool._blocks.size() - _node_pool._free.size();
        }
        
        [[nodiscard]] 
        uint32_t get_voxel_blocks_count() const
        {
            return _voxel_pool._blocks.size() - _voxel_pool._free.size();
        }

        void alloc_bulk(spatial_voxel* voxels, uint32_t count)
        {
            // would prefer to pre-alloc blocks, but the estimation can't be get right without proper 
            // voxel array spatial analysis and we ofc cant do that as it can potentially 
            // decrease performance in orders of magnitude.
            // Also, over estimating block counts even a bit will most of the time just
            // decrease performance (benchmarked on my own PC, specs info in documents), 
            // std::vector default allocation strategy seems to be the most consistent,
            // and generally the most optimal to rely on with bulk allocs for now
            for(uint32_t i = 0; i < count; ++i){
                vector_type pos{0,0,0}; 
                voxels[i].decode_position(&pos[0]);
                alloc(pos, voxels[i]._data);
            }
        }

        void alloc(vector_type& voxel_position, voxel_format& voxel)
        { 
            const bool overflow = 
            voxel_position[0] >= BOUNDS[0] || 
            voxel_position[1] >= BOUNDS[1] || 
            voxel_position[2] >= BOUNDS[2]; 

            if constexpr (DETAILS._discard_overflow){   
                if(overflow) { return; } 
            } else {
                assert(!overflow);
            }

            vector_type voxel_transformed; 

            // x2 even space
            voxel_transformed = voxel_position;
            voxel_transformed <<= 1; 
            
            node_format* node_ = &_root_node;
        
            // root node origin is always zero
            vector_type node_position{0,0,0}; 
            
            auto node_block_index = 0;

            ////////////////////////
            // TRAVERSE NODE TREE //
            ////////////////////////            
            // not tested with MSVC, with Clang on windows getting ~2x performance (my PC with 32^3 alloc benchmark, ~1.6m ns -> ~800k ns) 
            // improvement from succesful unrolling for this specific traversing loop, this indicates signifigant overhead is caused by 
            // branching and pipeline stalls
            #if defined(__clang__) || defined(__GNUC__)
            #pragma unroll
            #elif defined(_MSC_VER)
            #pragma loop(unroll)
            #endif
            for(int i = 0; i < MAX_DEPTH-2; ++i)
            {
                const auto depth_ = i;

                const auto node_extent = (component_type)(AXIS_WIDTH >> depth_);

                SVO_NODE_RELATION_MATH_IMPL();

                node_->_mask |= child_bit;

                node_block_index = node_->_block_index;
                auto& node_block_ = _node_pool._blocks[node_block_index];
        
                if(!exist)
                {   
                    node_format new_node_{};
                    new_node_._depth = depth_+1;
                    new_node_._mask  = 0;
                    new_node_._block_index = _node_pool.acquire_next_index();
                    node_block_[index] = new_node_;

                    _node_pool.alloc();

                    node_ = &_node_pool._blocks[node_block_index][index]; 
                } 
                else 
                {
                    node_ = &node_block_[index]; 
                }

                node_position = next_node_position;
            }

            ////////////////////////
            //    VOXEL OCTANT    //
            ////////////////////////
            {
                auto node_extent = (component_type)(4);

                SVO_NODE_RELATION_MATH_IMPL();

                node_->_mask |= child_bit;
                auto& node_block_ = _node_pool._blocks[node_->_block_index];

                if(!exist)
                {
                    node_format new_node_{};
                    new_node_._depth = MAX_DEPTH-1;
                    new_node_._mask  = 0;
                    // this time acquire voxel block index from voxel_pool
                    new_node_._block_index = _voxel_pool.acquire_next_index();
                    
                    node_block_[index] = new_node_;

                    // alloc new voxel block
                    _voxel_pool.alloc();

                    node_ = &node_block_[index]; 
                } 
                else 
                {
                    node_ = &node_block_[index];
                }
                node_position = next_node_position;
            }

            ////////////////////////
            //     ALLOC VOXEL    //
            ////////////////////////
            {
                auto node_extent = (component_type)(2);

                SVO_NODE_RELATION_MATH_IMPL();

                node_->_mask |= child_bit;

                auto& voxel_block_ = _voxel_pool._blocks[node_->_block_index];
                
                // allocate voxel
                voxel_block_[index] = voxel;
            }
        }

        voxel_format* get(const vector_type& voxel_position)
        {
            const bool overflow = 
            voxel_position[0] >= BOUNDS[0] || 
            voxel_position[1] >= BOUNDS[1] || 
            voxel_position[2] >= BOUNDS[2]; 

            if constexpr (DETAILS._discard_overflow){   
                if(overflow) { return nullptr; } 
            } else {
                assert(!overflow);
            }

            vector_type voxel_transformed; 

            // x2 even space
            voxel_transformed = voxel_position;
            voxel_transformed *= 2; 
            
            node_format* node_ = &_root_node;
        
            // root node origin is always zero
            vector_type node_position{0,0,0}; 

            ////////////////////////
            // TRAVERSE NODE TREE //
            ////////////////////////

            #if defined(__clang__) || defined(__GNUC__)
            #pragma unroll
            #elif defined(_MSC_VER)
            #pragma loop(unroll)
            #endif
            for(int i = 0; i < MAX_DEPTH-1; ++i) 
            {
                const auto depth_ = i;

                const auto node_extent = (component_type)(AXIS_WIDTH >> depth_);

                SVO_NODE_RELATION_MATH_IMPL();

                if(!exist){
                    return nullptr;
                }

                auto& node_block_ = _node_pool._blocks[node_->_block_index];
                node_ = &node_block_[index]; 
                node_position = next_node_position;
            };
            
            ////////////////////////
            //        VOXEL       //
            ////////////////////////

            auto node_extent = (component_type)(2);

            SVO_NODE_RELATION_MATH_IMPL();

            if(!exist)
            {
                return nullptr;
            }

            return &_voxel_pool._blocks[node_->_block_index][index];  
        }

        bool dealloc(const vector_type& voxel_position)
        {
            const bool overflow = 
            voxel_position[0] >= BOUNDS[0] || 
            voxel_position[1] >= BOUNDS[1] || 
            voxel_position[2] >= BOUNDS[2]; 

            if constexpr (DETAILS._discard_overflow){   
                if(overflow) { return false; } 
            } else {
                assert(!overflow);
            }

            uint8_t depth_ = 0;
            std::array<node_format*, MAX_DEPTH> path_{};
            std::array<uint8_t, MAX_DEPTH> child_bits{};
            auto voxel = 
                get_traced(
                    voxel_position, 
                    &path_[0], 
                    &child_bits[0], 
                    &depth_);
            
            if(!voxel) {
                return false;
            }
            
            auto& voxel_mask = path_[depth_]->_mask;
            voxel_mask &= ~child_bits[depth_];
            if(voxel_mask != 0){
                return true;
            }

            // # if we got here, we can dealloc the whole voxel block

            ///////////////////////////
            //  DEALLOC VOXEL BLOCK  //
            ///////////////////////////
            
            _voxel_pool.dealloc(path_[depth_]->_block_index);
            path_[depth_-1]->_mask &= ~child_bits[depth_-1];
            --depth_;
            
            ///////////////////////////
            //  DEALLOC NODE BLOCKS  //
            ///////////////////////////
            
            // traverse up, root node excluded
            for(; depth_ >= 1; --depth_){
                if(path_[depth_]->_mask != 0) { 
                    break; }
                _node_pool.dealloc(path_[depth_]->_block_index);
                path_[depth_-1]->_mask &= ~child_bits[depth_-1];
            }

            return true;
        }

        voxel_format* get_traced(const vector_type& voxel_position, node_format** node_path, uint8_t* child_bits, uint8_t* reached_depth)
        {
            const bool overflow = 
            voxel_position[0] >= BOUNDS[0] || 
            voxel_position[1] >= BOUNDS[1] || 
            voxel_position[2] >= BOUNDS[2]; 

            if constexpr (DETAILS._discard_overflow){   
                if(overflow) { return nullptr; } 
            } else {
                assert(!overflow);
            }

            vector_type voxel_transformed; 

            // x2 even space
            voxel_transformed = voxel_position;
            voxel_transformed *= 2; 
            
            node_format* node_ = &_root_node;
        
            // root node origin is always zero
            vector_type node_position{0,0,0}; 

            ////////////////////////
            // TRAVERSE NODE TREE //
            ////////////////////////

            #if defined(__clang__) || defined(__GNUC__)
            #pragma unroll
            #elif defined(_MSC_VER)
            #pragma loop(unroll)
            #endif
            for(int i = 0; i < MAX_DEPTH-1; ++i) 
            {
                const auto depth_ = i;

                const auto node_extent = (component_type)(AXIS_WIDTH >> depth_);

                SVO_NODE_RELATION_MATH_IMPL();

                if(!exist){
                    return nullptr;
                }

                ////////////////////////
                //        TRACE       //
                ////////////////////////

                child_bits[depth_] = child_bit;
                node_path [depth_] = node_;
                *reached_depth = depth_;

                auto& node_block_ = _node_pool._blocks[node_->_block_index];
                node_ = &node_block_[index]; 
                node_position = next_node_position;
            };

            ////////////////////////
            //        VOXEL       //
            ////////////////////////

            auto node_extent = (component_type)(2);

            SVO_NODE_RELATION_MATH_IMPL();

            if(!exist)
            {
                return nullptr;
            }

            ////////////////////////
            //        TRACE       //
            ////////////////////////

            constexpr auto DEPTH_END =  MAX_DEPTH-1;
            child_bits [DEPTH_END] = child_bit;
            node_path  [DEPTH_END] = node_;
            *reached_depth = DEPTH_END;

            return &_voxel_pool._blocks[node_->_block_index][index];  
        }
    };
}

