#!/bin/bash
if [ -d "build" ]; then
    read -p "Build directory already exists, would you like to delete it? (y/n): " choice
    case "$choice" in 
      y|Y ) rm -rf build; mkdir build;;
      n|N ) echo "Continuing with existing build directory...";;
      * ) echo "Invalid input, exiting..."; exit;;
    esac
else
    mkdir build
fi
cd build
cmake ..
make