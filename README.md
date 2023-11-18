# Moment and Weight Matrix Generation "In One Place"

## Introduction
This tool is designed for Simulated Method of Moments (SMM) estimation. The tool gives functionality to generate moments from the data, based on identification strategy specified in the configuration `.yml` file. The same function can then be used to generate the moments
from the simuldated data. The tool also generates the weight matrix based on influence functions of the *empirical* moments. 

Influence function generation follows [Erickson and Whited (2000) ](https://www.jstor.org/stable/10.1086/317670). See also Jay Kahn's useful exposition [here](https://j-kahn.com/files/influencefunctions.pdf). 

The project is currently under development and highly experimental.

## Features
- Handles panel data with multiple variables and time, cohort, and treatment variables.
- Computes mean, standard deviation, pairwise correlations, and autocorrelations.
- Supports gender and treatment group-based analysis.
- Outputs raw moments and influence functions for each moment. 

## Installation
Requires Pandas and YAML libraries. Install the required libraries using pip:

```bash
pip install pandas pyyaml
```


## License
This project is licensed under [MIT License](LICENSE).

