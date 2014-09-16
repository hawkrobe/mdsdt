mdsdt
=====

An R package for multi-dimensional signal detection theory, also known as GRT -- General Recognition Theory. 

We implement the recommendations from Silbert & Thomas (2013), providing up-to-date statistical methods to conduct GRT analyses. Following the result that failures of decisional and perceptual separability are generically confounded, we assume decisional separability and focus on testing perceptual separability and independence, providing the following options:

1. Full Gaussian model fitting with visualization and model comparison tools
2. Tests of marginal response invariance (MRI) to detect failures of perceptual separability .
3. Tests of sampling independence (SI) to detect failures of perceptual independence. 

Until the package is approved on CRAN, the easiest way to test it is to download the [source](https://iu.box.com/s/e5r2ljs70uwzl1r351dn), put it in the working directory, and install using the command

install.packages("mdsdt_0.1.tar.gz", repos = NULL, type = "source")

