# pdLISTA

Spatio-temporal cluster detection and classification through product density LISTA functions. This repository is based on `spatstat`, `splancs`, `stpp` and `KernSmooth` packages.

## Installation guide

The easiest way to install the development version of pdLISTA from github is using the devtools package which can be installed run the next command:
```
install.packages('devtools', dependencies=TRUE)
```
and thereafter run the commands:
```
require(devtools)
install_github('frajaroco/pdLISTA')
```
## References
- Siino, M., Rodríguez-Cortés, F. J., Mateu, J. and Adelfio, G. (2017). Testing for local structure in spatio-temporal point pattern data. **Submitted**.

## CiteBibtex
```
@misc{r16,
	author = {Francisco J. Rodr\'{i}guez-Cort\'{e}s},
	title = {pdLISTA: Second-order product density local indicator of spatio-temporal association function},
	year = {2016},
	note = {GitHub repository},
	url = {\url{https://github.com/frajaroco/pdLISTA}}
}
