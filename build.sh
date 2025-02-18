#!/bin/bash
mkdir -p build
pyinstaller --onefile TopMPI.py --distpath build/dist --workpath build/tmp --specpath build/spec
echo "Build complete! Executable is in build/dist/"

