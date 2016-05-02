# PMRF - Protein Markov random field
This is an evolutionary analysis tool using Markov random field.


## Installation

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

- PMRF is distributed under the MIT license. Please see the [file](LICENSE).
- PMRF includes the [libLBFGS] (currently version 1.10) library. [libLBFGS] is licensed under the MIT license. Please see the [file](thirdparty/lbfgs/lbfgs.h).
- PMRF includes the [Eigen] (currently version 3.2.8) library. [Eigen] is licensed under the MPL2 license. Please see the [file](thirdparty/eigen3/COPYING.README).


## Acknowledgements
This software includes and uses [libLBFGS] and [Eigen] libraries.


[PMRF repository]: https://github.com/jeongchans/pmrf/releases
[Eigen]: http://eigen.tuxfamily.org/
[libLBFGS]: http://www.chokkan.org/software/liblbfgs/
