# quasiNewtonMM
R package that implements the acceleration scheme for slowly-convergent MM algorithms (or any monotonous fixed-point iteration method) proposed by Agarwal and Xu (2020). It finds the fixed point of MM residual function using classical Broyden's quasi-Newton (BQN) method that makes secant approximations utilizing information from the MM algorithm map.

To install:
```{r}
# install.packages("devtools")
library(devtools)
devtools::install_github("medhaaga/quasiNewtonMM")
```
