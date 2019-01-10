## EnsDAM installation directory

To compile the library and the examples :

- create a 'make.macro' file corresponding to your compiler in the '../macro' directory.
  This is the Makefile configurable part, which specifies
  your compiler options and where to find the NetCDF library.

- edit the file 'Makefile.template' to include this 'make.macro' file (first line below the title)

- compile (library and examples) with:

```bash
make
make examples
```

To update the Makefile (if the source are modified) :

```bash
./mkmf -t Makefile.template -p ../lib/libensdam.a ../src/*
```