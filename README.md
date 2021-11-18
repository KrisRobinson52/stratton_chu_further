stratton-chu-cpp
==================

Implementation of Stratton-Chu integrals for electromagnetic beams propagation and reflection problem


# Building

To build the library you are required to build VTK from [VTK Repo](https://github.com/Kitware/VTK/tree/33519ffc861b46b3d780bc2ed47100af91b64dd3)
somewhere in your system.
To achieve this follow instructions README.md build step.
After that you can build the project with:
```bash
mkdir build && cd build
cmake .. -DVTK_DIR=<path to VTK build folder>
cmake --build . -j $(nproc)
```
You will find your executable in the `build/cpp/app/` with name `sc-runner`
