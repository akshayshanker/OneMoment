# Moment and Weight Matrix Generation "In One Place"

## Introduction
This tool is designed for Simulated Method of Moments (SMM) estimation. The tool gives functionality to generate moments from the data, based on identification strategy specified in the configuration `.yml` file. The exact same function can then be used to generate the moments
from the simuldated data. The tool also generates the weight matrix based on influence functions of the *empirical* moments. 

Influence function generation follows [Erickons and Whited (2000) ](https://www.jstor.org/stable/10.1086/317670). See also Jay Kahn's useful exposition [here](https://j-kahn.com/files/influencefunctions.pdf). 

The project is currently under development and highly experimental.

## Features
- Handles panel data with multiple variables and time, cohort, and treatment variables.
- Computes mean, standard deviation, pairwise correlations, and autocorrelations.
- Supports gender and treatment group-based analysis.
- Outputs raw moments and influence functions for each moment. 

## Usage

Read the config file which has identification strategy and data variable lists,  and read the data:

```
config = read_config('moments_LS.yml')
df_in = pd.read_csv('example_LS.csv')
``` 

Generate moments from the data. 

```
df, psi_df, raw_moments = generate_moments(df_in, config)
```

Compute the variance -covariance matrix of the moments. 

```
CovM = _compute_cov_matrix_with_nan(raw_moments[('Female', 'None')],
                                     psi_df[('Female', 'None')])
``` 


## License
This project is licensed under [MIT License](LICENSE).

