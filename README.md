<div align="center">
<h1><ins>GIIM</h1>
<h2>
    Generalized Incremental Image Mosaicking with a Coarse-to-Fine
    Framework via Graph Cuts
</h2>
</div>

## Usage

Download from:  
https://drive.google.com/drive/folders/1BeFTfLwQ7DjmYqA8Xhi7VUxj4ZgG00ig?usp=drive_link

### Building Tools

We have provided build tools for the Windows platform, including CMake configuration files and third-party libraries supplied by vcpkg.

### Pre-built Binaries

We also provide the pre-built binaries for the Windows platform.  
Run the following command to generate the seamline:

```bash
RS-Toolset_app_mosaicking.exe --input=<a folder for storing images> --output-mosaicking-vector=<a path ending with .shp> --low-overviews-trunc=3 --high-overviews-trunc=1
