# PMRF - Protein Markov random field
PMRF suite is a generic tool for Markov random field modeling and evolutionary analysis of protein sequences.


## Installation

### Requirements
PMRF suite does not have any dependencies other than the C++11 standard library.

### Download
Download the latest release from the [PMRF repository].

### Building the software
1. Uncompress the downloaded file in the directory where the PMRF software will be installed.

2. Compile the source code.
  ```
  $ make
  ```

3. The executable named as `pmrf` will be placed in the directory.


## Getting started

### Preparing input files
To generate an MRF model, a multiple alignment is required as an input. The multiple alignment can be formatted in FASTA format. As a default, PMRF uses a complete graph topology for the MRF model, but the customized MRF architecture can be used by determining the edge list.

PMRF suite considers the first aligned sequence in the multiple alignment as the reference sequence. Thus, the gapped positions in the first aligned sequences will ignored. In MRF model, the non-gapped positions will be considered as the nodes, and denoted by the position numbers starting with one.

The edge list consits of node pairs adjacent to each other. In edge list file, each line describes the neighboring nodes with their corresponding position numbers.

As an example, a [multiple alignment file](example/MYG_PHYCD.afa) and a [edge list file](example/MYG_PHYCD.edge) are provided in this package.

### Generating an MRF model from an MSA
For a given file, the `build` command carries out the parameterization for MRF weights, and writes the model into a file. By using the example files above, the MRF model can be built, and written to the model file, `MYG_PHYCD.mrf`, as below.

  ```
  $ pmrf build example/MYG_PHYCD.afa --edge example/MYG_PHYCD.edge -o MYG_PHYCD.mrf
  ```

If the option, `--edge example/MYG_PHYCD.edge`, is omitted, the output model archtecture will use the complete graph topology.

### Calculating the coevolution scores
The MRF model can be used for estimating evolutionary constraints. PMRF suite offers `stat` command to estimate the pairwise or positional coevolution scores.

The pairwise coevolution scores are calculated for edges of MRF model, by using `--mode pair` option.
  ```
  $ pmrf stat MYG_PHYCD.mrf --mode pair
  Pos1 Pos2      Score    Z-score
  ---- ---- ---------- ----------
     2    3   0.405721   1.644454
     2    4  -0.144241  -0.483282
     2   81  -0.152190  -0.514037
                :
  ```

The columns of output table represent the edge (Pos1-Pos2), the pairwise coevolution score (Score), and the Z-score transformation (Z-score), respectively.

The positional coevolution score is an alternative coevolution estimate determined for each node, while the pairwise coevolution score is detemined for each edge. As many biological informations are annotated per positions, this positional score would be useful. The positional coevolution scores are calculated by using `--mode pos` option.
  ```
  $ pmrf stat MYG_PHYCD.mrf --mode pos
   Pos      Score    Z-score
  ---- ---------- ----------
     1   0.205870   0.150840
     2   0.137707  -0.473971
     3   0.339084   1.371917
     4   0.200826   0.104600
     5   0.116322  -0.669993
     6   0.142767  -0.427584
     7   0.255363   0.604503
     8   0.198990   0.087774
               :
  ```

The columns of output table represent the node (Pos), the positional coevolution score (Score), and the Z-score transformation (Z-score), respectively.

### Estimating the statistical potential
Given a sequence, its statistical potential can be estimated with the likelihood from MRF model. PMRF suite provides the `infer` command to estimate the pseudo-likelihood values.

To estimate the statistical potential, a MRF model file and a sequence file formatted with FASTA are required. Please, note that the sequences should be aligned to the reference sequence of the MRF model, which implies that each of aligned sequences has the same length with the reference sequence length. Here, we prepared an [example file](example/MYG_SEQ.afa).
  ```
  $ pmrf infer MYG_PHYCD.mrf example/MYG_SEQ.afa
       -------- MRF -------- ------ Profile ------
     #      Score       Diff      Score       Diff Description
  ---- ---------- ---------- ---------- ---------- --------------------
     1     232.12       0.00     124.81       0.00 sp|P02185|MYG_PHYCD Myoglobin
     2     218.42     -13.69     103.16     -21.64 sp|P14399|MYG_MUSAN Myoglobin
     3     222.06     -10.05     105.93     -18.87 sp|P02206|MYG_HETPO Myoglobin
     4     235.05       2.94     112.74     -12.07 sp|P86874|MYG_DRONO Myoglobin
  ```

The columns of output table represent the sequence number (#), the pseudo-likelihood value (Score), its difference from the reference sequence (Diff), the profile-based likelihood value (Score), its difference from the reference sequence (Diff), and the description in the input FASTA file (Description), respectively.

### Retrieving the MRF parameters for further applications
Although the MRF model is writtein in a binary file, all the model paramters are easily accessible by using `show` command. To retrieve a node or edge parameter, use the `show` command with `--node <pos>` or `--edge <pos,pos>` option. For example, the node parameter corresponding to the sequence position `7` can be retrieved as below.
  ```
  $ pmrf show MYG_PHYCD.mrf --node 7
  A        0.2537
  C       -0.2907
  D       -0.3643
  E       -0.0548
  F       -0.1063
           :
  ```

The values represent the parameter values for the corresponding amino acids, respectively.

Similarly, the edge parameter interconnecting the sequence positions `7` and `134` can be retrieved as below.
  ```
  $ pmrf show MYG_PHYCD.mrf --edge 7,134
  7\134         A       C       D       E ...
  A        0.0441 -0.0588 -0.0069 -0.0150 ...
  C       -0.0926  0.0001  0.0140  0.0125 ...
  D       -0.1031  0.0062  0.0164  0.0152 ...
  E        0.0997  0.0640  0.0058  0.0032 ...
                        :
  ```

The values represent the parameter values for the corresponding amino acid pairs, respectively, and the row and column headers indicate the amino acids of the position `7` and `134`, respectively.

### Getting more help
PMRF suite offers vairous customization with command-line options. To find its capability and customize for your interests, see the help page.

The summarized command list is availalbe by `help` command.
  ```
  $ pmrf help
  Usage: pmrf <command> [<args>]
  
  The following commands are available for MRF modeling and varisous applications:
  
  generate MRF model
    build        Build MRF model
  
  examine evolutionary information
    infer        Estimate sequence distribution with MRF model
    stat         Estimate evolutionary constraints for MRF model
    show         Show MRF model parameters
  ```

To see the help page and the full list of options for a specific command, just run `pmrf <command>` with no arguments.

## Licenses
PMRF suite is distributed under the MIT license. Please see the [LICENSE](LICENSE) file for details.

Note that PMRF suite includes and uses third-party libraries licensed as below.

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
