
#include "rapid_svo.h"
#include <sstream>

#define ANKERL_NANOBENCH_IMPLEMENT
#include <nanobench.h>

#define TERMINAL_ANSI_RED(args) "\x1B[31m" << args << "\033[0m"
#define TERMINAL_ANSI_CYAN(args) "\x1B[36m" << args << "\033[0m"
#define TERMINAL_ANSI_MAGENTA(args) "\x1B[35m" << args << "\033[0m"
#define TERMINAL_ANSI_GREEN(args) "\x1B[32m" << args << "\033[0m"
#define TERMINAL_ANSI_YELLOW(args) "\x1B[33m" << args << "\033[0m"

static std::string add_thousand_separators(std::string value, char thousandSep = '.')
{
    int len = value.length();
    int dlen = 3;

    while (len > dlen) {
        value.insert(len - dlen, 1, thousandSep);
        dlen += 4;
        len += 1;
    }
    return value;
}

template<typename SVO_TREE_T>
static void svo_bench(std::stringstream& strbuf, std::string name, int extent, int max_epochs = 11, int max_iters = 50)
{
    using namespace ankerl::nanobench;
         
    auto svo_bench = ankerl::nanobench::Bench().minEpochIterations(50);
    svo_bench.output(nullptr);

    // VOXELS ALLOC / GET / DEALLOC
    {   
        using svo_tree = SVO_TREE_T;
        using svo_component_type = rapid_svo::morton_util<svo_tree::get_type()>::component_type;
        std::vector<typename svo_tree::spatial_voxel> voxels{};
        std::vector<typename svo_tree::spatial_voxel> voxels2{};
                
        auto max_depth_ = svo_tree::MAX_DEPTH;
                
        for(int x = 0; x < extent; ++x){
            for(int y = 0; y < (extent<2?1:extent/2); ++y){
                for(int z = 0; z < extent; ++z){
                    typename svo_tree::spatial_voxel voxel;
                    typename svo_tree::vector_type prr{
                        (svo_component_type)x,
                        (svo_component_type)y,
                        (svo_component_type)z};
                            
                    voxel.encode_position(&prr[0]); 
                    voxels2.push_back(voxel);
                }
            }
        }
        strbuf << "\n[" << TERMINAL_ANSI_CYAN(name) << "]";
        svo_tree tree{};
        tree.alloc_bulk(voxels2.data(), voxels2.size());
        strbuf << " \x1B[33m" << (tree.byte_size() / 1000.0) << " KB";
        strbuf << ", max_depth=" << svo_tree::MAX_DEPTH;
        strbuf << "\033[0m" << "\n";

        std::vector<typename svo_tree::vector_type> positions{};
        const auto count = extent*extent*extent; 
        for(int x = 0; x <extent; ++x){
            for(int y = 0; y < extent; ++y){
                for(int z = 0; z < extent; ++z){
                    positions.emplace_back(x,y,z);
                }
            }
        }

        for(int x = 0; x < extent; ++x){
            for(int y = 0; y < extent; ++y){
                for(int z = 0; z < extent; ++z){
                    typename svo_tree::spatial_voxel voxel;
                    typename svo_tree::vector_type pos_{
                        (svo_component_type)x,
                        (svo_component_type)y,
                        (svo_component_type)z};
                    voxel.encode_position(&pos_[0]); 
                    voxels.push_back(voxel);
                }
            }
        }

        using namespace ankerl::nanobench;
        svo_bench.epochs(max_epochs);
        svo_bench.minEpochIterations(max_iters);
        svo_bench.run("svo_alloc(N^3)", [&] {
            svo_tree tree2{};
            tree2.alloc_bulk(voxels.data(), voxels.size());
        });

        int found_voxels_count = 0; 
        for(int i = 0; i < count; ++i) {
            found_voxels_count += (tree.get(positions[i]) != nullptr);
        }

        svo_bench.run("svo_get(N^3)", [&] {
            for(int i = 0; i < count; ++i) {
                doNotOptimizeAway(tree.get(positions[i]));
            }
        });

        svo_bench.epochs(1);
        svo_bench.minEpochIterations(1);
        svo_bench.run("svo_dealloc(N^3)", [&] {
            for(int i = 0; i < count; ++i) {
                doNotOptimizeAway(tree.dealloc(positions[i]));
            }
        });
    }

    const auto count = extent*extent*extent; 
    auto& results = svo_bench.results();
    for(int i = 0; i < results.size(); ++i)
    {
        auto r = results[i];
        auto time_per_op_ns = r.median(ankerl::nanobench::Result::Measure::elapsed) * 1e9;
        double op_per_s = 1e9 / time_per_op_ns; 
        uint32_t voxel_count = op_per_s * count;
        strbuf << std::fixed << add_thousand_separators(std::to_string(voxel_count)) << " voxels/s\t| ";
        strbuf << std::fixed << std::setprecision(2) << time_per_op_ns << " ns/op\t| ";
        strbuf << std::fixed << std::setprecision(2) << (time_per_op_ns/count) << " ns/voxel\t| ";
        strbuf << std::fixed << op_per_s << " op/s\t| ";
        strbuf << r.config().mBenchmarkName << " as " << count << " voxels/op";
        strbuf << "\n";
    }
}

int main()
{
    using namespace ankerl::nanobench;
    
    printf("\nRAPID_SVO running benchmarks... this will take a moment\n");

    std::stringstream
    all_results{};

    constexpr rapid_svo::details_info details_16b_16pow3{
        ._discard_overflow = true,
        ._limit_max_bounds = { 16,16,16 }};

    constexpr rapid_svo::details_info details_16b_32pow3{
        ._discard_overflow = true,
        ._limit_max_bounds = { 32,32,32 }};

    constexpr rapid_svo::details_info details_32b_32pow3{
        ._discard_overflow = true,
        ._limit_max_bounds = { 32,32,32 }};
           
    constexpr rapid_svo::details_info details_32b_64pow3{
        ._discard_overflow = true,
        ._limit_max_bounds = { 64,64,64 }};

    constexpr rapid_svo::details_info details_32b_128pow3{
        ._discard_overflow = true,
        ._limit_max_bounds = { 128,128,128 }};

    constexpr rapid_svo::details_info details_32b_256pow3{
        ._discard_overflow = true,
        ._limit_max_bounds = { 256,256,256 }};

    constexpr rapid_svo::details_info details_32b_512pow3{
        ._discard_overflow = true,
        ._limit_max_bounds = { 512,512,512 }};

    constexpr rapid_svo::details_info details_32b_1024pow3{
        ._discard_overflow = true,
        ._limit_max_bounds = { 1024,1024,1024 }};

    constexpr rapid_svo::details_info details_32b_64pow3_full_space{};

    constexpr rapid_svo::details_info details_32b_128pow3_full_space{};

    constexpr rapid_svo::details_info details_32b_256pow3_full_space{};

    std::stringstream strbuf{};

    //all_results << svo_bench<svo::tree<svo::morton_32b, svo::basic_voxel_format, details_32b_1024pow3>>
    //("32b_space__svo_bench(128^3)__64bit_voxels", 128, 1,1);
            
    svo_bench<rapid_svo::tree<rapid_svo::morton_16b, rapid_svo::basic_voxel_format, details_16b_16pow3>>
    (strbuf, "16b_space__svo_bench(16^3)__64bit_voxels", 16);

    svo_bench<rapid_svo::tree<rapid_svo::morton_16b, rapid_svo::basic_voxel_format, details_16b_32pow3>>
    (strbuf, "16b_space__svo_bench(32^3)__64bit_voxels", 32);

    svo_bench<rapid_svo::tree<rapid_svo::morton_32b, rapid_svo::basic_voxel_format, details_32b_32pow3>>
    (strbuf, "32b_space__svo_bench(32^3)__64bit_voxels", 32);

    svo_bench<rapid_svo::tree<rapid_svo::morton_32b, rapid_svo::basic_voxel_format, details_32b_64pow3>>
    (strbuf, "32b_space__svo_bench(64^3)__64bit_voxels", 64, 6, 10);

    svo_bench<rapid_svo::tree<rapid_svo::morton_32b, rapid_svo::basic_voxel_format, details_32b_128pow3>>
    (strbuf, "32b_space__svo_bench(128^3)__64bit_voxels", 128, 1, 1);
            
    std::cout << strbuf.str() << std::flush;

    return 1;
}