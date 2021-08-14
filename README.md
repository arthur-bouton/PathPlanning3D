Tested on Ubuntu 18.04 and requires OpenGL and cmake 3.10+.

## Installation

if you don't have OpenGL installed:

```sudo apt install xorg-dev libglu1-mesa-dev freeglut3-dev```

Then build the application:

```mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
./main```
