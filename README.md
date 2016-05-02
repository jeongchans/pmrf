# PMRF - Protein Markov random field
This is an evolutionary analysis tool using Markov random field.


## Installation

### Requirements
PMRF does not have any dependencies other than the C++11 standard library.

### Download
Download the latest release from the [PMRF repository].

### Building the software
1. Uncompress the downloaded file in the directory where the PMRF software will be installed.

2. Compile the source code.
  ``
  $ make
  ``
3. The executable will be placed in the directory.


## Usage
TBA

### Help messages
Show the general help page.
  ``
  $ pmrf --help
  ``
Show the help page for a specific command.
  ``
  $ pmrf <COMMAND> --help
  ``


## Examples

### Generating MRF model from an MSA
TBA

### Generating MRF model by specifying the network architecture
`pmrf-build` can specify the network architecture of MRF model with `--edge` option.

### Calculating the coevolution scores for an MRF model
TBA

### Calculating the pseudolikelihood of aligned sequeces
TBA


## Licenses
PMRF is distributed under the MIT license. Please see the [LICENSE](LICENSE) file for details.

Note that PMRF includes and uses third-party libraries licensed as below.

  - libLBFGS
    Files: thirdparty/lbfgs/*
    License: MIT license
    Please see the [thirdparty/lbfgs/lbfgs.h](thirdparty/lbfgs/lbfgs.h) file for details.

  - Eigen
    Files: thirdparty/eigen/*
    License: MPL2 license
    Please see the [thirdparty/eigen3/COPYING.README](thirdparty/eigen3/COPYING.README) file for details.


## Acknowledgements
This software includes and uses [libLBFGS] (currently version 1.10) and [Eigen] libraries (currently version 3.2.8).


[PMRF repository]: https://github.com/jeongchans/pmrf/releases
[libLBFGS]: http://www.chokkan.org/software/liblbfgs/
[Eigen]: http://eigen.tuxfamily.org/
