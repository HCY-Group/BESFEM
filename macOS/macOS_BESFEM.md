## macOS Installation

BESFEM can be compiled on macOS using Homebrew. The macOS Makefile uses Homebrew to locate MFEM, Open MPI, Hypre, SuiteSparse, and libTIFF.

### 1. Install Homebrew

If Homebrew is not already installed, install it from the official Homebrew website.

Verify the installation with:

```bash
brew --version
```

### 2. Install Required Dependencies

Install the packages required by BESFEM:

```bash
brew install git
brew install open-mpi
brew install mfem
brew install hypre
brew install suite-sparse
brew install libtiff
```

The Makefile uses `brew --prefix` to determine where these packages are installed, so no hard-coded `/opt/homebrew` or `/usr/local` paths are required.

You can verify the installations with:

```bash
brew --prefix mfem
brew --prefix open-mpi
brew --prefix hypre
brew --prefix suite-sparse
brew --prefix libtiff
```

For example, on an Apple Silicon Mac these will typically be located under:

```text
/opt/homebrew/opt/
```

### 3. Clone BESFEM

Clone the BESFEM repository:

```bash
git clone https://github.com/HCY-Group/BESFEM.git
cd BESFEM
```

### 4. Compile BESFEM

From the root BESFEM directory, run:

```bash
make clean
make -f macOS/Makefile.mac
```

The resulting executable will be created at:

```text
bin/battery_simulation
```
