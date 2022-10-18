## CrownDensity
CrownDensity simulates growth of one Scots pine tree growing in a gap of a homogeneous forest. The density of the surrounding forest is input from file stemsha.fun. CrownDensity uses the code of project LignumForest for tree growth.

## Compilation
This project must reside under *lignum-core* directory. That means first clone lignum-core repository and then in lignum-core clone [CrownDensity] (https://github.com/lignumsystem/CrownDensity.git). You will of course need to have also LignumForest in lignum-core.

To compile CrownDensity (and lignum-core) type:

    cd CrownDensity
    qmake  CrownDensity.pro
    make

To compile with optimization on (faster, no debug) type:

    qmake  "CONFIG+=release" CrownDensity.pro
    make

To remove all compilation work type `make distclean`.

## Running the program

Command line options and their short documentation can be obtained by running the program
without any parameters: <CODE> ./crowndens </CODE>
You will will see the output from usage(): \sa GrowthLoop::usage()


## Documentation

The introductionary presentation is in [GENERAL_DESCRIPTION.md](GENERAL_DESCRIPTION.md).


The Reference Guide for the CrownDensity will be based on comments and other information
available in the software. Extraction of the comments, rendition of the software content and 
architecture, generation of the structure of the document and formatting the document to html 
and LaTeX will be done by `doxygen`. To generate the documentation run `doxygen` in CrownDensity directory:
    
    doxygen Doxyfile 2> errors.txt
     
Doxyfile is the configuration file for `doxygen`. The documentation will appear in DoxygenDoc directory. 
Errors and warnings will appear in *errors.txt*. To see html version of the document type (on macOS):

    open DoxygenDoc/html/index.html
    
To generate LaTeX version go to latex subdirectory and use make:

    cd DoxygenDoc/latex
    make all
    
The result will be *refman.pdf* that can be opened with a pdf reader.

To use Doxyfile the following three programs are needed:

  + doxygen: generate the document 
  + dot: used by `doxygen` to generate graphs for class hierarchies and function calls.
  + doxywizard: GUI to browse, edit and optionally run Doxyfile. 
    
On macOS these are easiest to install with MacPorts (or some other software package system). 
