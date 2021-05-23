# Progetto-CUDA

Problema di linking (credo)

in: mcc solveC.tm solveC.cpp -o solveC

out: /usr/lib/gcc/x86_64-linux-gnu/6/../../../x86_64-linux-gnu/Scrt1.o: In function `_start':
(.text+0x20): undefined reference to `main'
/tmp/ccpTwiVy.o: In function `__static_initialization_and_destruction_0(int, int)':
solveC.cpp:(.text+0x74d): undefined reference to `std::ios_base::Init::Init()'
solveC.cpp:(.text+0x762): undefined reference to `std::ios_base::Init::~Init()'
/tmp/cc8fAmCL.o: In function `_tr0':
solveC.tm.c:(.text+0x3d): undefined reference to `solveC'
collect2: error: ld returned 1 exit status

dim: 10000
mathematica error:
No more memory available.
Mathematica kernel has shut down.
Try quitting other applications and then retry.
