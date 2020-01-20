
## SIRUS: Stable and Interpretable RUle Set
C. Benard (Safran Tech & Sorbonne University)

### Introduction
SIRUS (Stable and Interpretable RUle Set) is a regression and classification algorithm based on random forests, which takes the form of a short list of rules. 
SIRUS combines the simplicity of decision trees with the predictivity of random forests for problems with low order interactions. 
The core aggregation principle of random forests is kept, but instead of aggregating predictions, SIRUS selects the most frequent nodes of the forest to form a stable rule ensemble model.

This R package is a fork from the  project ranger (https://github.com/imbs-hl/ranger).

### Installation
To install sirus, download sirus-0.2.1.tar.gz and just run:
	install.packages("sirus-0.2.1.tar.gz", repos = NULL)

R version >= 3.1 is required. If you compile yourself, the new RTools toolchain is required.

Pre-compiled binaries are available for windows (sirus_0.2.1.zip), mac (sirus_0.2.1.tgz) and linux (sirus_0.2.1_R_x86_64-redhat-linux-gnu.tar.gz).


### Usage
For usage of the package see the Examples section. 

If you find any bug, or have any suggestions for improvement, please let us know !

### References
* Benard C., Biau G., da Veiga S., Scornet E. (2019). SIRUS: making random forests interpretable. <arXiv:1908.06852>. https://arxiv.org/abs/1908.06852
* Breiman, L. (2001). Random forests. Mach Learn, 45:5-32. https://doi.org/10.1023/A:1010933404324
* Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. https://doi.org/10.18637/jss.v077.i01

