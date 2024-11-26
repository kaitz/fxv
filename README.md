# fxv

fxv is a file archiver and compressor development program based on a paq8 architecture. Programm uses config files for compression, detection, decoding and encoding.
Compression models and decoding models are saved into the final archive. Decompressor uses them when data is extracted.
Main compression routine is stored uncompressed at the start of an archive.

fxv uses process virtual machine, which compiles [FX language](help/FXlanguage.md) source code to bytecode at runtime and executes it. 
Command line option -j allows the use of x86 as JIT target.

## Usage example
__fxv.exe  -1 -j data__

Compress using x86 JIT mode (j) config file is __config.pxv__ by defaul and input file is __data__.


# Contents
[Overview](#overview)

[Command line options](help/CommandLineOptions.md)

[FX language](help/FXlanguage.md)

[Config files overview](help/ConfigFiles.md)

[Example Config files](help/ExampleConfigs.md)

[More Examples and Tuning](help/MoreExamples.md)

[Components](help/components.md)

[Forum](#forum)

[History](#history)

[Testing results](#testing-results)

# Overview
![image](https://github.com/user-attachments/assets/e4a728db-22fe-4b1f-a3ff-de3356475caf)

# History
See paq8pxv [paq8pxv]( https://github.com/kaitz/paq8pxv )

# Testing results
[Google spreadsheet file](https://docs.google.com/spreadsheets/d/1vWunOvKT5hMmy6BQB6MIofJ6AL2GJqXxXOGbC_W_-WQ/edit?usp=sharing)
