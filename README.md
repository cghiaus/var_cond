# Variable conductivity

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/var_cond/HEAD)

## Introduction

In general, [thermal conductivity](https://en.m.wikipedia.org/wiki/Thermal_conductivity_and_resistivity) varies with temperature, $\lambda = \lambda(\theta)$. The variation may be modelled by a linear or a quadratic function.

When the variation is linear, the model is

$$
\lambda = \lambda_0 [1 + \beta (\theta - \theta_0)]
$$

or

$$\lambda = \lambda_0 (1 + \beta \theta)$$

if the change of variable $(\theta - \theta_0) \rightarrow \theta$ is done.

When the variation is quadratic, the model is

$$
\lambda = \lambda_0 [1 + \beta (\theta - \theta_0)^2]
$$

or

$$\lambda = \lambda_0 (1 + \beta \theta^2)$$

if the change of variable $(\theta - \theta_0) \rightarrow \theta$ is done.

The coefficient $\beta$ is called _temperature coefficient of thermal conductivity_ (Cengel & Ghajar, 2020).

## Research questions

We try to answer the reseach questions :

1. Does the conductivity variation with temperature influences significantly the results as compared with the use of a constant value?

2. If the answer is yes, how can we integrate effectivelly the this variation in models.


## References

__Theory__
- Cengel, Y. A., & Ghajar, A. J. (2015). Heat and Mass Transfer: Fundamentals and Applications, 5th Edition, McGraw-Hill Education. New York. ISBN 978-0-07-339818-1

- [Cengel Y. & Ghajar A. (2020)](https://www.studocu.com/in/document/priyadarshini-engineering-college/english/htchapter-02-xyz/42524065). Chapter 2 Heat condution equatiion, in Heat and Mass Transfer: Fundamentals and Applications, 6th Edition, McGraw-Hill Education, ISBN10: 0073398195 | ISBN13: 9780073398198

- [Cengel Y. & Ghajar A. (2020a)](https://www.studocu.com/row/document/celal-bayar-universitesi/engineering-mechanics/heat-chap02-094-this-is-summaries/11179160). Chapter 2 Heat condution equatiion, in Heat and Mass Transfer: Fundamentals and Applications, 6th Edition, Solution manual, McGraw-Hill Education

- [Berardi U. et al.(2018)](https://doi.org/10.3390/en11040872). On the Effects of Variation of Thermal Conductivity in
Buildings in the Italian Construction Sector, Energies, 11(4), 872

- [Wang, Y., Zhang, S., Wang, D., & Liu, Y. (2023)](https://doi.org/10.1016/j.enbenv.2022.02.008). Experimental study on the influence of temperature and humidity on the thermal conductivity of building insulation materials. Energy and Built Environment, 4(4), 386-398

- [Higgis, B.G. (2020)](https://doi.org/10.13140/RG.2.2.17178.36805) Heat Conduction with Variable Thermal Conductivity, University of California, Davis

__SI Units__
- [BIPM (2019)](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf/2d2b50bf-f2b4-9661-f402-5f9d66e4b507?version=1.11&t=1671101192839&download=true). The International System of Units (SI), 9th edition, licence CC-BY-3.0

- [Gőbel, E., Mills, I., Wallard,  A. (2006)](https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-concise-EN.pdf/2fda4656-e236-0fcb-3867-36ca74eea4e3). A concise summary of the International System of Units, the SI

- [Thomson, A., Taylor, B. N. (2008)](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication811e2008.pdf). Guide for the use of the international System of units (NIST Special Publication 811․ 2008 Edition). National Institute of Standards and Technology, US Government Printing Office.

__Scientific computing__

- [Wilson, G., et al.(2014)](https://doi.org/10.1371/journal.pbio.1001745). Best practices for scientific computing. PLoS biology, 12(1), e1001745.

- [Wilson, G., et al. (2017)](https://doi.org/10.1371/journal.pcbi.1005510). Good enough practices in scientific computing. PLoS computational biology, 13(6), e1005510.

- [Noble, W.S. (2009)](https://doi.org/10.1371/journal.pcbi.1000424) A Quick Guide to Organizing Computational Biology Projects. PLoS Comput Biol 5(7): e1000424.

- [Balaban, G. et al (2021)](https://doi.org/10.1371/journal.pcbi.1008549) Ten simple rules for quick and dirty scientific programming, PLoS Comput Biol 17(3): e1008549

- [Cookiecutter Data Science](https://drivendata.github.io/cookiecutter-data-science/), GitHub

- [Organization and Packaging of Python Projects](https://rabernat.github.io/research_computing/organization-and-packaging-of-python-projects.html)