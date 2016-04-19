# PMRF - Protein Markov random field
This is an evolutionary analysis tool using Markov random field.


## Installation

### Requirements
The [Eigen] library is essential for building and running this software. We used the [Eigen], version 3.2.8.

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


## Acknowledgements
This software includes and uses [libLBFGS] (http://www.chokkan.org/software/liblbfgs/).


[PMRF repository]: https://github.com/jeongchans/pmrf/releases
[Eigen]: http://eigen.tuxfamily.org/
[libLBFGS]: http://www.chokkan.org/software/liblbfgs/
