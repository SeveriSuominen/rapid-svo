```text
 ██████╗  █████╗ ██████╗ ██╗██████╗     ███████╗██╗   ██╗ ██████╗        
 ██╔══██╗██╔══██╗██╔══██╗██║██╔══██╗    ██╔════╝██║   ██║██╔═══██╗       
 ██████╔╝███████║██████╔╝██║██║  ██║    ███████╗██║   ██║██║   ██║       
 ██╔══██╗██╔══██║██╔═══╝ ██║██║  ██║    ╚════██║╚██╗ ██╔╝██║   ██║       
 ██║  ██║██║  ██║██║     ██║██████╔╝    ███████║ ╚████╔╝ ╚██████╔╝       
 ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝     ╚══════╝  ╚═══╝   ╚═════╝
```
## Description 

Single-header fast dynamic sparse voxel octree implementation using integer coordinates and Morton codes. This implementation is the result of prototyping solutions for lightning-fast voxel physics simulations on CPU in chunked procedural game worlds. As a technical challenge, I have been micro-optimizing this to the best of my abilities, but I am confident that further optimizations can be achieved. 

I get considerable faster results by compiling with Clang compared to MSVC, I have not yet determined which optimizations are taking place to account for the difference. Benchmarks seems to be very consistent with both Clang and MSVC.

## Building

CMake can be used to build the project, all depedencies are fetched automatically.

Build: ``` cmake . -B build ```

Run: ``` cmake --build build --config Release ```

## Specs
My personal setup used with the benchmarks: **AMD Ryzen 7 1800X Eight-Core, 32GB DDR4, Windows 10**

## Benchmarks 

### Clang 19.1.6
```text
[16b_space__svo_bench(16^3)__64bit_voxels] 18.888 KB, max_depth=4
58.277.015 voxels/s     | 70285.00 ns/op        | 17.16 ns/voxel        | 14227.79 op/s | svo_alloc(N^3) as 4096 voxels/op
176.495.371 voxels/s    | 23207.41 ns/op        | 5.67 ns/voxel         | 43089.69 op/s | svo_get(N^3) as 4096 voxels/op
116.970.950 voxels/s    | 35017.24 ns/op        | 8.55 ns/voxel         | 28557.36 op/s | svo_dealloc(N^3) as 4096 voxels/op

[16b_space__svo_bench(32^3)__64bit_voxels] 149.96 KB, max_depth=5
52.990.425 voxels/s     | 618375.86 ns/op       | 18.87 ns/voxel        | 1617.14 op/s  | svo_alloc(N^3) as 32768 voxels/op
102.509.747 voxels/s    | 319657.41 ns/op       | 9.76 ns/voxel         | 3128.35 op/s  | svo_get(N^3) as 32768 voxels/op
103.141.328 voxels/s    | 317700.00 ns/op       | 9.70 ns/voxel         | 3147.62 op/s  | svo_dealloc(N^3) as 32768 voxels/op

[32b_space__svo_bench(32^3)__64bit_voxels] 149.96 KB, max_depth=5
47.332.906 voxels/s     | 692287.93 ns/op       | 21.13 ns/voxel        | 1444.49 op/s  | svo_alloc(N^3) as 32768 voxels/op
113.550.730 voxels/s    | 288575.86 ns/op       | 8.81 ns/voxel         | 3465.29 op/s  | svo_get(N^3) as 32768 voxels/op
107.789.473 voxels/s    | 304000.00 ns/op       | 9.28 ns/voxel         | 3289.47 op/s  | svo_dealloc(N^3) as 32768 voxels/op

[32b_space__svo_bench(64^3)__64bit_voxels] 1198.54 KB, max_depth=6
34.022.648 voxels/s     | 7704985.00 ns/op      | 29.39 ns/voxel        | 129.79 op/s   | svo_alloc(N^3) as 262144 voxels/op
98.313.496 voxels/s     | 2666409.09 ns/op      | 10.17 ns/voxel        | 375.04 op/s   | svo_get(N^3) as 262144 voxels/op
50.198.962 voxels/s     | 5222100.00 ns/op      | 19.92 ns/voxel        | 191.49 op/s   | svo_dealloc(N^3) as 262144 voxels/op

[32b_space__svo_bench(128^3)__64bit_voxels] 9587.14 KB, max_depth=7
29.914.485 voxels/s     | 70104900.00 ns/op     | 33.43 ns/voxel        | 14.26 op/s    | svo_alloc(N^3) as 2097152 voxels/op
89.402.195 voxels/s     | 23457500.00 ns/op     | 11.19 ns/voxel        | 42.63 op/s    | svo_get(N^3) as 2097152 voxels/op
51.225.512 voxels/s     | 40939600.00 ns/op     | 19.52 ns/voxel        | 24.43 op/s    | svo_dealloc(N^3) as 2097152 voxels/op
```

### MSVC 1937 (VS 17.7)
```text
[16b_space__svo_bench(16^3)__64bit_voxels] 18.888 KB, max_depth=4
22.732.612 voxels/s     | 180181.67 ns/op       | 43.99 ns/voxel        | 5549.95 op/s  | svo_alloc(N^3) as 4096 voxels/op
42.565.425 voxels/s     | 96228.33 ns/op        | 23.49 ns/voxel        | 10391.95 op/s | svo_get(N^3) as 4096 voxels/op
84.403.885 voxels/s     | 48528.57 ns/op        | 11.85 ns/voxel        | 20606.42 op/s | svo_dealloc(N^3) as 4096 voxels/op

[16b_space__svo_bench(32^3)__64bit_voxels] 149.96 KB, max_depth=5
20.362.022 voxels/s     | 1609270.37 ns/op      | 49.11 ns/voxel        | 621.40 op/s   | svo_alloc(N^3) as 32768 voxels/op
34.347.387 voxels/s     | 954017.24 ns/op       | 29.11 ns/voxel        | 1048.20 op/s  | svo_get(N^3) as 32768 voxels/op
31.100.987 voxels/s     | 1053600.00 ns/op      | 32.15 ns/voxel        | 949.13 op/s   | svo_dealloc(N^3) as 32768 voxels/op

[32b_space__svo_bench(32^3)__64bit_voxels] 149.96 KB, max_depth=5
19.614.832 voxels/s     | 1670572.55 ns/op      | 50.98 ns/voxel        | 598.60 op/s   | svo_alloc(N^3) as 32768 voxels/op
34.128.711 voxels/s     | 960130.00 ns/op       | 29.30 ns/voxel        | 1041.53 op/s  | svo_get(N^3) as 32768 voxels/op
29.925.114 voxels/s     | 1095000.00 ns/op      | 33.42 ns/voxel        | 913.24 op/s   | svo_dealloc(N^3) as 32768 voxels/op

[32b_space__svo_bench(64^3)__64bit_voxels] 1198.54 KB, max_depth=6
15.255.098 voxels/s     | 17184025.00 ns/op     | 65.55 ns/voxel        | 58.19 op/s    | svo_alloc(N^3) as 262144 voxels/op
28.903.534 voxels/s     | 9069617.27 ns/op      | 34.60 ns/voxel        | 110.26 op/s   | svo_get(N^3) as 262144 voxels/op
26.166.515 voxels/s     | 10018300.00 ns/op     | 38.22 ns/voxel        | 99.82 op/s    | svo_dealloc(N^3) as 262144 voxels/op

[32b_space__svo_bench(128^3)__64bit_voxels] 9587.14 KB, max_depth=7
12.638.516 voxels/s     | 165933400.00 ns/op    | 79.12 ns/voxel        | 6.03 op/s     | svo_alloc(N^3) as 2097152 voxels/op
24.787.302 voxels/s     | 84605900.00 ns/op     | 40.34 ns/voxel        | 11.82 op/s    | svo_get(N^3) as 2097152 voxels/op
22.847.704 voxels/s     | 91788300.00 ns/op     | 43.77 ns/voxel        | 10.89 op/s    | svo_dealloc(N^3) as 2097152 voxels/op
```
