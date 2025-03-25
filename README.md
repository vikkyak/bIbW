# bIbW

**bIbW**: *Between- and Within-sample Normalization for Single-cell RNA-seq Data Using Partial Least Squares and Adaptive Fuzzy Weights*

This R package provides a normalization method for single-cell RNA-seq data that reduces variation across cells and conditions using partial least squares (PLS) regression and adaptive fuzzy weighting.

---

## 📦 Installation

You can install the development version of `bIbW` from GitHub:

```r
# install.packages("devtools")
devtools::install_github("vikkyak/bIbW")


🛠️ Usage

library(bIbW)

# Simulate data
set.seed(123)
G <- 1000; n <- 300
mu <- rgamma(G, shape = 2, rate = 2)
NB_cell <- function(j) rnbinom(G, size = 0.1, mu = mu)
SimulatedData <- sapply(1:n, NB_cell)

# Run normalization
Clusters <- rep(1:2, each = 150)
Result <- bIbWNoRm(Data = SimulatedData, Clusters = Clusters, Method = "KernSmooth")

📖 Documentation

Run ?bIbWNoRm in R after installing to access the full documentation.

 Author
Vikas Singh
Email: vikkysingh07@gmail.com

License
This package is distributed under the MIT license.


---

### 🔹 3. Add, commit, and push

```bash
git add README.md
git commit -m "Add README.md to describe the package"
git push origin main











