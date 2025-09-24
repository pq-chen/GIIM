<div align="center">
<h1><ins>GIIM</h1>
<h2>
    Generalized Incremental Image Mosaicking with a Coarse-to-Fine
    Framework via Graph Cuts
</h2>
</div>

## Usage

Google Drive Link:  
https://drive.google.com/drive/folders/1BeFTfLwQ7DjmYqA8Xhi7VUxj4ZgG00ig?usp=drive_link

### Building Tools

We provide build tools for the Windows platform, including CMake configuration files and third-party libraries managed through vcpkg.

### Pre-built Binaries

Pre-built binaries for the Windows platform are also available, together with four GaoFen-2 images from the paper as test data. To generate the seamline, run the following command:

```bash
RS-Toolset_app_mosaicking.exe \
  --input=<Guangzhou-GF2 directory> \
  --output-mosaicking-vector=<a path ending with .shp> \
  --low-overviews-trunc=3 \
  --high-overviews-trunc=1
