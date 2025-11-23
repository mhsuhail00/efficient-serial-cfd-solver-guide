# Serial CFD Solver Modernization (v0 – v5)
### Recommendations on writing Serial Computational Fluid Dynamic Code

[![Paper](https://img.shields.io/badge/Research%20Paper-FMFP%202025-blue)](https://doi.org/your-doi-here)
[![Language](https://img.shields.io/badge/C++-17-00599C?logo=cplusplus)](https://isocpp.org/)
[![FORTRAN](https://img.shields.io/badge/FORTRAN-77-734F96)](https://fortran-lang.org/)
[![CFD](https://img.shields.io/badge/CFD-Computational%20Fluid%20Dynamics-green)]()
[![Memory](https://img.shields.io/badge/Optimization-Memory%20Layout-red)]()
[![Blitz++](https://img.shields.io/badge/Library-Blitz%2B%2B-yellow)]()
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE) 

## Overview

This repository contains the complete implementation of six progressively optimized CFD solver versions (v0–v5) developed for the research paper:

**"Recommendations on Writing Serial Computational Fluid Dynamic Code"**  
*Faran Imam, M. Suhail, Abdus Samad, and Syed Fahad Anwer*  
*Proceedings of the 12th International and 52nd National Conference on Fluid Mechanics and Fluid Power (FMFP)*  
*December 19-21, 2025, Nirma University, Gujarat, India*<br>
**DOI**: [https://doi.org/your-doi-here](https://doi.org/your-doi-here)

This project modernizes a legacy FORTRAN 77 finite-difference CFD solver into a sequence of optimized C++ (C++17) implementations. Six versions (v0–v5) progressively improve performance by redesigning memory layout, data access, and abstraction layers.

The repository includes:

- **solver_0 — FORTRAN 77 Legacy Solver**  
- **solver_1 — Direct C++ translation**  
- **solver_2 — Loop-reordered for row-major locality**  
- **solver_3 — Fully flattened memory model (1D arrays + Index Increment)**  
- **solver_4 — Blitz++ Array Allocation**  
- **solver_5 — Blitz++ Array allocation + manual pointer traversal (Reduce call overheads)**  

All versions solve the same 2D incompressible flow benchmark and are intended for **performance comparison, teaching, and reproducible research**.


## Requirements (Linux)

The following tools and libraries are required to build and benchmark the solvers:

- **gfortran**: GNU Fortran compiler for v0 (FORTRAN 77)
- **g++**: GNU C++ compiler with C++17 support
- **Blitz++**: Scientific array library (required for v4 and v5)
- **Python 3**: For performance visualization (matplotlib, numpy)
- **make**: Build automation tool
- **perf**: Linux performance monitoring tool (optional, for profiling)

### Installing Dependencies

**Ubuntu/Debian:**
```bash
# Install compilers and build tools
sudo apt-get update
sudo apt-get install gfortran g++ make

# Install Blitz++ library
sudo apt-get install libblitz0-dev

# Install Python visualization dependencies
sudo apt-get install python3 python3-pip
pip3 install matplotlib numpy

# Install perf (for performance profiling)
sudo apt-get install linux-tools-common linux-tools-generic
```

## Execution Guide
### Quick Start

To run the complete benchmark suite (compile, execute, and visualize):

```bash
make
```

This single command executes the entire pipeline automatically.

### Execution Workflow

#### Step 1: Automated Benchmarking (`make run`)

The `execute.sh` script performs the following operations for each solver version (v0–v5):

1. **Execution with Profiling**: Runs each solver for `MAXSTEP` iterations (default: 5×10⁶)
   - Collects performance metrics using `perf stat`
   - Measures: execution time, cache misses, instruction count, branch operations
   
2. **Data Collection**: Saves profiling results to `./temp_dir/`
   - Creates text files: `sol0_perf.txt`, `sol1_perf.txt`, ..., `sol5_perf.txt`
   - Each file contains: CPU time, cache statistics, instruction counts

**Note**: The default iteration count is 5×10⁶ (5 million), which may take several hours to days depending on your system. To reduce execution time for testing:

```cpp
// Edit in each solver file before running
MAXSTEP = 1000;  // Reduce to 1000 iterations for quick testing
```

#### Step 2: Visualization (`make plot`)

The `plot_metrics.py` script:

1. **Reads** performance data from `./temp_dir/*.txt`
2. **Generates** comparison plots:
   - Execution time vs. solver version
   - Cache miss rates (L1/L2/L3)
   - Instruction overhead
3. **Saves** figures to `./plots/` directory as PNG/PDF files

### Individual Makefile Targets

```bash
# Run only the benchmark (without plotting)
make run

# Generate plots from existing data
make plot

# Clean temporary performance data
make clean

# Clean all generated files (data + plots)
make clean_all

# Make execute.sh executable (done automatically by 'make run')
make prepare
```

### Manual Execution (Without Makefile)

If you prefer to run solvers individually:

```bash
# Compile and run FORTRAN baseline
gfortran solver_0.for -o solver_0
./solver_0

# Compile and run C++ version 3 (flattened memory)
g++ solver_3.cpp -o solver_3
./solver_3

# Compile and run version 5 (Blitz++ with pointers)
g++ solver_5.cpp -o solver_5 -lblitz
./solver_5
```

## Results Summary

Performance over 5×10⁶ iterations on Intel i9-12900K:

|           Solver Versions         |    CPU Time    |  Speedup  |    Cache Efficiency   |
|-----------------------------------|----------------|-----------|-----------------------|
| solver_0 (FORTRAN)                |   8.866 days   |   1.00×   |      Baseline         |
| solver_1 (C++ Exact Conversion)   |   26.606 days  |   0.33×   |   Poor (fragmented)   |
| solver_2 (Reordered Nested Loops) |   10.113 days  |   0.88×   |      Improved         |
| solver_3 (Flattened Array)        |   5.444 days   |   1.63×   |      Excellent        |
| solver_4 (Blitz Arrays)           |   54.549 days  |   0.16×   | Poor (Context Switch) |
| solver_5 (Blitz + Pointer Access) | **4.626 days** | **1.92×** |      **Optimal**      |

**Recommendation**: Version 5 (Blitz++ allocation + pointer traversal) provides the best balance of performance, maintainability, and memory safety.

## Authors

**Mohammad Suhail**  
Department of Computer Engineering  
Zakir Husain College of Engineering and Technology (ZHCET)  
Aligarh Muslim University, Aligarh 202002, India  
mhsuhail0@outlook.com

**Faran Imam**  
Department of Computer Engineering  
Zakir Husain College of Engineering and Technology (ZHCET)  
Aligarh Muslim University, Aligarh 202002, India  
faranimam4@gmail.com

## Acknowledgments

The authors extend their gratitude to **Prof. Syed Fahad Anwer** and **Dr. Abdus Samad** for their mentorship, encouragement, and constructive guidance throughout this research.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.



