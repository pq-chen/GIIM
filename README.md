<div align="center">
<h1><ins>GIIM</h1>
<h2>
    Generalized Incremental Image Mosaicking with a Coarse-to-Fine
    Framework via Graph Cuts
</h2>
</div>

## Installation

All files referenced below can be downloaded from the following Google Drive link:  
https://drive.google.com/drive/folders/1BeFTfLwQ7DjmYqA8Xhi7VUxj4ZgG00ig?usp=drive_link

### Build from Source

Build tools for the Windows platform are provided, including CMake configuration files and third-party libraries managed via vcpkg.

1. Download **LCMake** and **vcpkg-export-20250521-153211.zip**, and extract them.  
2. Specify the toolchain file for cross-compiling located in  
   `vcpkg-export-20250521-153211\scripts\buildsystems\vcpkg.cmake`.  
3. Set `LCMAKE_DIR` to `LCMake\cmake`.  
4. Configure and generate the project.  

### Pre-built Binaries

Pre-built binaries for the Windows platform are available in **GIIM_release.zip**.

## Usage

Four GaoFen-2 images used in the paper are provided as test data in **Guangzhou-GF2.zip**.  
To generate the seamline, run the following command:

```bash
RS-Toolset_app_mosaicking.exe \
  --input=<Guangzhou-GF2 directory> \
  --output-mosaicking-vector=<a path ending with .shp> \
  --low-overviews-trunc=3 \
  --high-overviews-trunc=1
