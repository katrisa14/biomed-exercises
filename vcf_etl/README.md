# VCF Transformation

VCF files come in a variety of formats.  We often have to take vendor formats and transform them in small and large ways to create a common VCF format for our system.  An example VCF file is provided in the "/data" directory.  To normalize this file 4 changes need to be done to it
1. Remove any records that don't have a FILTER of "PASS" or "LowGQX"
1. The INFO column should only contain allele frequency, i.e. "AF={actual value}"
1. The FORMAT column should only include the following values "GT:AD:DP"
1. The final file should only have variant rows with an alternate allele.

This wrapper project uses [PDM](https://pdm.fming.dev/latest/) Please reference that site for installation.  

Once properly installed and coded, users should be able to run
1. "pdm run process normalize {in-file} {out-file}" to execute the scripts
1. "pdm run test" to run all unit tests

Please create a branch off of this repo and submit your code as a PR.  If you're unable to create a github account, you may also zip up the project and send it to your point of contact.