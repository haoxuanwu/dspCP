## Test environments
- R-hub windows-x86_64-devel (r-devel)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
On windows-x86_64-devel (r-devel), fedora-clang-devel (r-devel)

❯ checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Haoxuan Wu <hw399@cornell.edu>'
  
  New submission
  
Explanation: This is new package we are uploading to CRAN.

❯ checking for detritus in the temp directory ... NOTE
  Found the following files/directories:
    'lastMiKTeXException'
    
Explanation: As noted in R-hub issue #503, this could be due to a bug/crash in MiKTeX and can likely be ignored.
    
0 errors ✔ | 0 warnings ✔ | 2 notes ✖
