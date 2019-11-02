<img src="README_files/img/dpe.png" align="left" height="80" vspace="8" hspace="18">

**psilence** is a package for the R statistical language, that implements differentially private algorithms for releasing statistics while safeguarding the privacy of individuals in the raw data.

It provides the underlying algorithms for *PSI (&Psi;): a Private data Sharing Interface*, a system created to enable researchers in the social sciences and other fields to share and explore privacy-sensitive datasets with the strong privacy protections of differential privacy ([Project](http://psiprivacy.org/about)|[Paper](https://arxiv.org/abs/1609.04340)).  However, it can be used by itself in any application needing code for differentially private mechanisms.

In order to use the `Snapping Mechanism`, you will need to have [python3](https://www.python.org/downloads/) installed, as well as have all libraries included in [snapping_requirements.txt](snapping_requirements.txt). The included [Dockerfile](Dockerfile) can be used to create a Docker container in which to run **PSIlence**, or you can use it as a blueprint to ensure you have the software you need to run **PSIlence** outside of a container.

PSI and the PSIlence package have been developed by the [Privacy Tools Project](http://privacytools.seas.harvard.edu), an interdisciplinary research group at Harvard, MIT and Boston University.
