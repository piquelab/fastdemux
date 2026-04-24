# fastdemux htslib build/install fix

## Summary
This repository had a build/install issue caused by hardcoded htslib paths in CMake.

Observed errors included:
- compile-time header failure: `htslib/sam.h: No such file or directory`
- install/runtime failure to resolve `libhts`

This fix makes htslib detection dynamic and keeps installation local to this project by default.

## Root cause
The previous CMake config used fixed paths for htslib include/lib directories. On systems where `samtools` is provided by a module (for example Longleaf), htslib can be installed under a different layout, such as:

- `<samtools_prefix>/htslib/include/htslib/sam.h`
- `<samtools_prefix>/htslib/lib/libhts.so`

When paths did not match, compile failed. Even when compile worked in-session, installed binaries could still fail if runtime library lookup was not configured.

## What was changed
The following changes were made in CMake:

1. Dynamic htslib discovery
- Uses `find_program(samtools)` to infer the active module prefix.
- Uses `find_path` for `htslib/sam.h`.
- Uses `find_library` for `libhts`.
- Supports both common layouts:
  - `include`, `lib`, `lib64`
  - `htslib/include`, `htslib/lib`

2. Correct target include/link behavior
- Adds discovered `HTSLIB_INCLUDE_DIR` to targets.
- Links targets to discovered `HTSLIB_LIBRARY`.

3. Reliable runtime linking after install
- Sets per-target `BUILD_RPATH` and `INSTALL_RPATH` to htslib library directory.
- Keeps `CMAKE_INSTALL_RPATH_USE_LINK_PATH` enabled.

4. Safe default install location
- If install prefix is not explicitly set, defaults to `<repo>/install`.
- This avoids changing global or user-wide defaults.

## Files changed
- `CMakeLists.txt`

## Build and install (recommended)
From the fastdemux repo root:

```bash
module load samtools

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=<install-dir>

cmake --build build -j
cmake --install build
```

## Verify installation
Check that binaries exist:

```bash
ls -l <install-dir>/bin/fastdemux
ls -l <install-dir>/bin/fastaseq
```

Check dynamic library resolution:

```bash
ldd <install-dir>/bin/fastdemux | grep -i hts
ldd <install-dir>/bin/fastaseq | grep -i hts
```

Expected result: `libhts.so` resolves to your samtools/htslib module path.

## Run fastdemux after install
Use the installed binary directly:

```bash
<install-dir>/bin/fastdemux \
  -t 2 \
  <bam> \
  <vcf.gz> \
  <barcodes.tsv.gz> \
  <output_prefix>
```

Required runtime inputs:
- BAM file
- indexed VCF (`.vcf.gz` + index)
- barcode file
- output prefix

## Notes on defaults and safety
This fix is project-scoped and does not modify system defaults.
- No global config files were changed.
- No shell startup files were changed.
- Install target is local to this project unless overridden.
