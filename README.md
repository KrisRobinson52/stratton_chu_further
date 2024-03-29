stratton-chu-cpp
==================

Implementation of Stratton-Chu integrals for electromagnetic beams propagation
and reflection problem.


# Building

To build the library you are required to build VTK from [VTK Repo](https://github.com/Kitware/VTK/tree/33519ffc861b46b3d780bc2ed47100af91b64dd3)
somewhere in your system.
To achieve this follow instructions README.md build step.

Do not forget about Git submodules.
Clone the repo providing `--recurse-submodules` flag or if you've simply
cloned the repo in nonrecursive run these commands:
```bash
git submodule init
git submodule update --init --recursive
```

Download fftw3 library from [FFTW3 Link](http://fftw.org/fftw3_doc/Installation-and-Customization.html)
and install it by provided instructions

After that you can build the project with:
```bash
mkdir build && cd build
cmake .. -DVTK_DIR=<path to VTK build folder>
cmake --build . --parallel $(nproc)
```
You will find your executable in the `build/cpp/app/` with name `sc-runner`
