## CrownDensity

CrownDensity simulates growth of one Scots  pine tree growing in a gap
of  a homogeneous  forest. The  density of  the surrounding  forest is
input from file *stemsha.fun*.   CrownDensity uses components from the
LignumForest       for      tree       growth.      Most       notably
LignumForest::ScotsPineTree,     LignumForest::ScotsPineSegment    and
LignumForest::SetScotsPineSegmentLength.

CrownDensity project  must reside under *lignum-core*  directory. That
means clone first [lignum-core](https://github.com/lignumsystem/lignum-core.git)
repository  and  then  in  *lignum-core*  clone [CrownDensity](https://github.com/lignumsystem/CrownDensity.git)
and [LignumForest](https://github.com/lignumsystem/LignumForest.git) 
repositories.

## Compilation

It  seems  that  [CMake](https://cmake.org) is  becoming  increasingly
popular  (perhaps   de  facto  standard)  build   system  for  various
programming languages including C/C++. This is true also for Qt system
that seems to be giving up its `qmake` tool for CMake (Qt6 onwards). 

The Lignum  core system as  well as two recent  projects, CrownDensity
and LignumForest, have now CMake project files to organise and compile
Lignum  software.  Qt  `qmake`  can be  still used  but  it is  highly
recommended to start to use CMake.

One of the main benefits of CMake is its ability to use separate build
trees, i.e.   the software compilation  takes under a  single separate
build directory located  outside the source trees.  The second benefit
is the  ability to generate  build systems for several  IDEs including
Xcode.

### CrownDensity compilation with CMake

To create CMake build system for  Unix Makefile build system for debug
(debug information)  first create  the build  tree directory  and then
with `cmake` the Unix Makefile build system itself. Under CrownDensity
type:

    mkdir debug
    cd  debug
    cmake .. -DCMAKE_BUILD_TYPE=Debug
	make install 

Optionally explicitely set C++ compiler (macOS example):

    cmake ..  -DCMAKE_BUILD_TYPE=Debug  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS=-stdlib=libc++
    make install

For Unix Makefile build system for Release (optimised, no debug information) type:

    mkdir release
    cd release
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make install

In all cases `make install` will move `crowndens` to CrownDensity directory
where there are two example  shell scripts to run the program:
	
	run-crowndens-basic-model.sh	
	run-crowndens.sh

For Xcode IDE create the Xcode project file:

    mkdir xcode
    cd xcode
    cmake .. -G Xcode

Open  Xcode  IDE  from  Terminal  (or  open  the  Xcode  project  file
`crowndens.xcodeproj` from XCode):
     
	 open crowndens.xcodeproj

Build the `crowndens` Product in  Xcode for debugging.  It will appear
in *xcode/Debug*  directory:

	Xcode -> Product (in the menu bar) -> Build For -> Running/Testing/Profiling

See  also that: 

	Xcode -> Product (in the menu bar) -> Scheme 

is set  to `crowndens` to allow Run: 

	Xcode -> Product (in the menu bar) -> Run
	
to debug the program.

Copy necessary \*.fun  function files and \*.txt parameter files to
*xcode/Debug*  where   the  `crowndens`  is  located   in  this  case.
Otherwise  hard coded  files names  in the  program are  not found  by
`crowndens`. You can also copy `crowndens` to CrownDensity project
directory instead and load the binary to Xcode from there. 

Set command  line parameters for  `crowndens` in 

	Xcode -> Product (in the menu  bar) -> Scheme ->  Edit Scheme -> Arguments.

Divide the command line into practical parts for debugging from `Arguments -> '+'`.

**NB1:** To let Unix Makefile build system recognise file dependencies
correctly  (for example  changes made  in  c++adt in  the Lignum  core
system) clean  first the build  tree from previous software  build. In
the build directory type:

	make clean
	make install
	
By  default  CMake  only  tracks   changes  made  in  current  project
(ie. CrownDensity). Xcode IDE itself tracks file dependencies.

**NB2:** To  remove all  compilation work just  remove the  build tree
directory (i.e. *debug*, *release* or *xcode*) and recreate  the build
tree directory.

**NB3:**   CMake  projects   are   configured  with   *CMakeLists.txt*
files.  The  best  way  to   learn  CMake  is  by  studying  examples.
lignum-core and CrownDensity provide  CMakeLists.txt file examples how
to create libraries,  find and integrate external  libraries and setup
the final product with its dependenices.


### CrownDensity compilation with qmake

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
You will will see the output from CrownDensity::Usage(): \sa CrownDensity::Usage()


## Documentation

The Reference Guide for the CrownDensity will be based on comments and other information
available in the software. Extraction of the comments, rendition of the software content and 
architecture, generation of the structure of the document and formatting the document to html 
and LaTeX will be done by `doxygen`. To generate the documentation run `doxygen` in CrownDensity directory:
    
    doxygen Doxyfile 2> errors.txt
     
Doxyfile is the  configuration file for `doxygen` and  `2>` (zsh, bash
shells) redirects erros to  *errors.txt* file.  The documentation will
appear in DoxygenDoc  directory.  To see html version  of the document
type (on macOS):

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

### LaTeX typesetting
Documentation uses LaTeX typesetting for equations. To generate LaTeX equations `amsmath` package must be part of LaTeX installation.
On macOS this can be done for example with MacPorts TeXLive full installation:

   	 sudo port install texlive +full.

Doxyfile has the necessary LaTeX `amsmath` configuration. 
