#!/bin/sh

mkdir -p build/
emcmake cmake -B build -S .
cmake --build build --target rsiscool

cp build/rsiscool.js dist/
cp build/rsiscool.d.ts dist/
cp build/rsiscool.wasm dist/

mkdir -p build-workers/
emcmake cmake -B build-workers -S .
cmake --build build-workers --target rsiscool_workers

cp build-workers/rsiscool_workers.js dist/
cp build-workers/rsiscool_workers.d.ts dist/
ln -f build/rsiscool.wasm dist/rsiscool_workers.wasm
