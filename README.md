# PMRF - Protein Markov random field
This is an evolutionary analysis tool using Markov random field.


## Installation

### Prerequisites
The [Blitz++] (http://blitz.sourceforge.net) library is essential for running this software.

### Download
Download the latest release from [the git repository] (https://github.com/jeongchans/pmrf/releases).

### Building the software
1. Uncompress the downloaded file in the directory where the PMRF software will be installed.

2. Compile the source code.
  ```sh
  $ make
  ```

3. The executable will be placed in the directory.


## Usage
TBA

### Help messages
Show the general help page.
  ```sh
  $ pmrf --help
  ```

Show the help page for a specific command.
  ```sh
  $ pmrf <COMMAND> --help
  ```


## Examples

### Generating MRF model from an MSA
TBA

### Generating MRF model by specifying the network architecture
`pmrf-build` can specify the network architecture of MRF model with `--edge` option.

### Calculating the pseudolikelihood for aligned sequeces
TBA


## Acknowledgements
This software includes and uses [libLBFGS] (http://www.chokkan.org/software/liblbfgs/).
