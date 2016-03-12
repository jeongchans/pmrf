*** Settings ***
Library         pmrf_lib.py


*** Test Cases ***
print help messages
    [Template]  Run
    ${EMPTY}    ${HELP MESSAGE}
    --help      ${HELP MESSAGE}
    build -h    ${BUILD HELP}
    infer -h    ${INFER HELP}
    stat -h     ${STAT HELP}

print version
    Run     --version       pmrf version 0.2.0

build MRF model
    Run     build ${AFA FILE} --edge ${EDGE FILE} -o model.mrf    ${EMPTY}
    Should exist    model.mrf
    Files should be same    model.mrf   ${MRF FILE}
    Remove file     model.mrf

build MRF model with A3M MSA
    Run     build ${A3M FILE} --msa a3m --edge ${EDGE FILE} -o model.mrf      ${EMPTY}
    Files should be same    model.mrf   ${MRF FILE}
    Remove file     model.mrf

calculate likelihood
    Run     infer ${MRF FILE} ${AFA FILE}   [.]*

calculate positional coevolution scores
    Run     stat ${MRF FILE} --mode pos     [.]*

calculate pairwise coevolution scores
    Run     stat ${MRF FILE} --mode pair     [.]*


*** Variables ***
${HELP MESSAGE} =       Usage: pmrf [.]*
${BUILD HELP} =         Usage: pmrf build [.]*
${INFER HELP} =         Usage: pmrf infer [.]*
${STAT HELP} =          Usage: pmrf stat [.]*

${A3M FILE} =           ../../results/20160310/MYG_PHYCD/MYG_PHYCD.a3m
${AFA FILE} =           ../../results/20160310/MYG_PHYCD/MYG_PHYCD.afa
${EDGE FILE} =          ../../results/20160310/MYG_PHYCD/MYG_PHYCD.edge
${MRF FILE} =           ../../results/20160310/MYG_PHYCD/MYG_PHYCD.mrf


*** Keywords ***
Run
    [Arguments]     ${opts}     ${expected}
    ${output} =     Run PMRF    ${opts}
    Should match regexp     ${output}   ${expected}

File should exist
    [Arguments]     ${filename}
    Should exist    ${filename}

Files should be same
    [Arguments]     ${filename1}    ${filename2}
    ${content1} =       Read file   ${filename1}
    ${content2} =       Read file   ${filename2}
    Should be equal     ${content1}     ${content2}
