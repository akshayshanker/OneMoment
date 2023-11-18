# Moment and Weight Matrix Generation "In One Place"

## Introduction
This tool is designed for Statistical Method of Moments (SMM) estimation. The tool gives functionality to generate moments from the data, based on identification strategy specified in the configuration `.yml` file. The exact same function can then be used to generate the moments
from the simuldated data. The tool also generates the weight matrix based on influence functions of the *empirical* moments. 

Influence function generation follows Toni Whited's [paper](https://www.jstor.org/stable/10.1086/317670). See also Jay Kahn's useful exposition [here](https://j-kahn.com/files/influencefunctions.pdf). 

The project is currently under development and highly experimental.

## Features
- Generates statistical moments for SMM estimation from complex datasets.
- Handles data with multiple variables, including time, cohort, and treatment variables.
- Computes mean, standard deviation, pairwise correlations, and autocorrelations.
- Supports gender and treatment group-based analysis.
- Outputs raw moments and influence functions for specified demographic groups.

## Installation
To use this tool, you'll need Python installed on your system, along with the Pandas and YAML libraries. Install the required libraries using pip:

```bash
pip install pandas pyyaml
```


## License
This project is licensed under [MIT License](LICENSE).

