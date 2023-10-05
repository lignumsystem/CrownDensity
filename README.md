# CrownDensity

CrownDensity simulates growth of one Scots  pine tree growing in a gap
of  a homogeneous  forest. The  density of  the surrounding  forest is
input from file *stemsha.fun*.   CrownDensity uses components from the
LignumForest for  tree growth. Most notably
   
   + LignumForest::ScotsPineTree
   + LignumForest::ScotsPineSegment
   + LignumForest::SetScotsPineSegmentLength

## CrownDensity dependencies

CrownDensity depends on Lignum core system  and [LignumForest](https://github.com/lignumsystem/LignumForest.git). 
CrownDensity and LignumForest projects  must reside under *lignum-core*  project. 
Clone first [lignum-core](https://github.com/lignumsystem/lignum-core.git)
repository  and  then  in  *lignum-core*  clone 
[CrownDensity](https://github.com/lignumsystem/CrownDensity.git)
and [LignumForest](https://github.com/lignumsystem/LignumForest.git) 
repositories.

## CMake build system

It  seems  that  [CMake](https://cmake.org) is  becoming  increasingly
popular (perhaps  de facto standard) cross-platform  build system for
various programming languages  including C/C++. This is  true also for
Qt system that seems  to be giving up its `qmake`  tool for CMake (Qt6
onwards).

The Lignum  core system as  well as two recent  projects, CrownDensity
and LignumForest, have now CMake project files to organise and compile
Lignum  software.  Qt  `qmake`  can be  still used  but  it is  highly
recommended to start to use CMake.

One of the main benefits of CMake is its ability to use separate build
trees from  source file  trees, i.e.   the software  compilation takes
under a  single separate  build directory  located outside  the source
trees.  The second benefit is  the ability to generate build processes
for  traditional  Unix  Makefile  systems   as  well  as  for  several
Integrated   Development  Environments   (IDE)  including   Xcode  and
Microsoft  Visual  Studio  from   the  same  set  of  *CMakeLists.txt*
configuration files.

## CrownDensity: CMake for macOS and Unix/Linux Makefile build system

To create Makefile build system with CMake first create the
build tree  directory and  then with `cmake`  the Unix  Makefile build
system itself. To build the Lignum core system:

	cd lignum-core
	mkdir build
	cd build 
	cmake .. 
	make install
	
See also *lignum-core* [README](https://github.com/lignumsystem/lignum-core/blob/master/README.md).

To create CrownDensity Makefile build system for debug and compile `crowndens` binary 
type:

    cd CrownDensity
    mkdir debug
    cd  debug
    cmake .. -DCMAKE_BUILD_TYPE=Debug
    make install 

For CrownDensity Makefile build system for Release (optimised, no debug information) type:

    cd CrownDensity
    mkdir release
    cd release
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make install

In both cases `make install` will move `crowndens` to CrownDensity directory
where there are two example  shell scripts to run the program:
	
    run-crowndens-basic-model.sh	
    run-crowndens.sh

Command line options and their  short documentation can be obtained by
running `./crowndens`  without any  command line parameters.  See also
CrownDensity::Usage().

>[!IMPORTANT]
>It is important to type `make install` to also move `crowndens` to
>directory above to be used by the scripts to run simulatations.
>Typing just `make` the `crowndens` program remains in the compilation directory.

>[!IMPORTANT]
>To let Unix Makefile build system keep up with file dependencies
>correctly  (for example  changes made  in  c++adt in  the Lignum  core
>system) clean  first the build  tree from previous software  build.

To recompile `crowndens` type:

	make clean
	make install
	
By default Unix Makefile build system tracks only changes made
in the current CrownDensity project.

>[!IMPORTANT]
>To remove all CMake  configurations and compilation work just
>remove the build  tree directory (i.e. *debug*,  *release* or *xcode*)
>and recreate the build tree directory.

>[!NOTE]
>CMake  projects   are   configured  with   *CMakeLists.txt*
>files. For  this CMake  has an  extensive set  of CMake  variables and
>built-in functions that can be set in CMakeLists.txt files or given in
>command line.

The best way to  learn CMake is by  studying examples.
lignum-core and CrownDensity provide  CMakeLists.txt file examples how
to create libraries, find and integrate external libraries (Qt, HDF5),
create and use external binaries (`l2c` to compile L-system files) and
setup the final product with its dependenices.

## CrownDensity: CMake for Xcode

For Xcode IDE create the Xcode project file:

    mkdir xcode
    cd xcode
    cmake .. -G Xcode

Open  Xcode  IDE  from  Terminal. Alternatively open  the  Xcode  project  file
`crowndens.xcodeproj` from XCode:
     
	 open crowndens.xcodeproj

Build the `crowndens` Product in  Xcode for debugging.  It will appear
in *xcode/Debug*  directory:

	Xcode -> Product (in the menu bar) -> Build For -> Running/Testing/Profiling

See  also that: 

	Xcode -> Product (in the menu bar) -> Scheme 

is set  to `crowndens` to allow Run: 

	Xcode -> Product (in the menu bar) -> Run
	
to debug the program. Xcode IDE itself tracks file dependencies.

Copy necessary \*.fun  function files and \*.txt parameter files to
*xcode/Debug*  where   the  `crowndens`  is  located   in  this  case.
Otherwise  hard coded  files names  in the  program are  not found  by
`crowndens`. You can also copy `crowndens` to CrownDensity project
directory instead and load the binary to Xcode from there. 

Set command  line parameters for  `crowndens` in Xcode:

	Xcode -> Product (in the menu  bar) -> Scheme ->  Edit Scheme -> Arguments.

Divide the command line into practical parts for debugging from `Arguments -> '+'`.

## CMake for CrownDensity dependency graph

CMake allows to generate `graphviz` output file to show all library and executable dependencies of the project.
Then with `dot` create image file with desired file format. For example in the *release* directory type:
	
	mkdir graphviz
	cmake ..   --graphviz=graphviz/CrownDensity.dot
	dot -Tpdf -Kneato -Goverlap=prism  graphviz/CrownDensity.dot  -o  CrownDensity.pdf
	
The output file *CrownDensity.pdf* contains the visual presentation of the target dependenices including
external binaries and required link libraries. The option `-T` understands many well known image file formats.


## CrownDensity compilation with qmake

To compile CrownDensity (and lignum-core) type:

    cd CrownDensity
    qmake  CrownDensity.pro
    make

To compile with optimization on (faster, no debug) type:

    qmake  "CONFIG+=release" CrownDensity.pro
    make

To remove all compilation work type `make distclean`.


## CrownDensity Documentation

The Reference Guide for the CrownDensity will be based on comments and
other  information  available  in  the  software.  Extraction  of  the
comments,  rendition   of  the  software  content   and  architecture,
generation  of  the  structure  of the  document  and  formatting  the
document to html and LaTeX will  be done by `doxygen`. To generate the
documentation run `doxygen` in CrownDensity directory:
    
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
  + graphviz: generate dependency graphs for function calls and class dependencies
  + dot: generate figures from graphviz files
  + doxywizard: GUI to browse, edit and optionally run Doxyfile 
    
On macOS these are easiest to install with MacPorts (or some other software package system). 

### LaTeX typesetting
Documentation uses LaTeX typesetting for equations. To generate LaTeX equations `amsmath` package must be part of LaTeX installation.
On macOS this can be done for example with MacPorts TeXLive full installation:

	sudo port install texlive +full.

Doxyfile for CrownDensity has the necessary LaTeX `amsmath` configuration. 

