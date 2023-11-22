"""
This module provides functionality to compute statistical moments and their 
variance-covariance matrix for Simulated Method of Moments (SMM) estimation.

It includes functions to calculate mean, standard deviation,
pairwise correlations, and autocorrelation for specified
variables in a DataFrame by time and age group.

This is an experimental module and is not yet fully tested.
Author: Akshay Shanker, University of New South Wales, Sydney, Australia
Date: 18 November, 2023
Email: akshay.shanker@me.com 

See main code for usage example.
"""

import numpy as np
import pandas as pd
import yaml

def _nan_cov(x, y):
    """
    Compute the covariance between two arrays while excluding NaN values.

    Parameters:
    ----------
    x : 1D array
        First array for covariance calculation.
    y : 1D array
        Second array for covariance calculation.

    Returns:
    -------
    float
        Covariance between x and y, excluding NaN values. Returns NaN if
        there are no valid (non-NaN) pairs.
    """
    valid_indices = ~np.isnan(x) & ~np.isnan(y)
    if np.any(valid_indices):
        return np.cov(x[valid_indices], y[valid_indices])[0, 1]
    else:
        return np.nan
        
def _compute_cov_matrix_with_nan(dict_moments, dict_influence_functions):
    """
    Compute the variance-covariance matrix of moments using influence functions.
    Observations with NaN values are excluded from the calculations.

    Parameters:
    ----------
    dict_moments : dict
        Dictionary with real-valued moments as keys.
    dict_influence_functions : dict
        Dictionary where keys are moment labels and values are influence 
        functions (1D arrays).

    Returns:
    -------
    numpy.ndarray
        The variance-covariance matrix of the moments.

    Notes:
    -----
    Moments and influence functions must have the same labels for correct matching.
    """

    # Extract the order of moments from the dictionary keys
    moment_order = list(dict_moments.keys())

    # Gather corresponding influence functions in the order of moments
    influence_matrix = np.array([
        dict_influence_functions[moment] for moment in moment_order
        if moment in dict_influence_functions
    ])

    # Initialize an empty matrix for the covariance calculations
    num_moments = len(moment_order)
    cov_matrix = np.full((num_moments, num_moments), np.nan)

    # Compute the covariance matrix, handling NaNs
    for i in range(num_moments):
        for j in range(i, num_moments):
            cov_matrix[i, j] = cov_matrix[j, i] = \
                _nan_cov(influence_matrix[i], influence_matrix[j])

    return cov_matrix


def _computeAC(df, variable, current_time, 
                                        previous_time, time_var, member_id_var):
    """
    Compute the autocorrelation for a variable across specified time periods
    in a DataFrame. It includes NA for individuals not present in both periods.
    The length of psi matches the number of individuals in the current period.

    Parameters:
    df : pandas.DataFrame
        The DataFrame containing the data.
    variable : str
        The variable for which autocorrelation is computed.
    current_time : str/int
        The current time period.
    previous_time : str/int
        The previous time period.
    time_var : str
        Column name representing time in the DataFrame.
    member_id_var : str
        Column name representing member IDs in the DataFrame.

    Returns:
    corr: float
        The autocorrelation coefficient.
    psi: numpy.ndarray
        The autocorrelation influence function.
    """
    # Select data for current and previous time periods
    current_data = df[df[time_var] == current_time][[member_id_var, variable]]
    previous_data = df[df[time_var] == previous_time][[member_id_var, variable]]

    # Merge data on member_id for aligned individuals across years
    merged_data = pd.merge(current_data, previous_data, on=member_id_var, 
                           suffixes=('_current', '_previous'), how='left')

    # Calculate autocorrelation for valid pairs
    valid_indices = merged_data[f'{variable}_previous'].notna()
    if valid_indices.any():
        corr = np.corrcoef(
            merged_data.loc[valid_indices, f'{variable}_current'],
            merged_data.loc[valid_indices, f'{variable}_previous'])[0, 1]

        # Calculate Psi for current members
        psi_values = (
            (merged_data[f'{variable}_current'] - 
             merged_data[f'{variable}_current'].mean()) *
            (merged_data[f'{variable}_previous'].fillna(0) - 
             merged_data.loc[valid_indices, f'{variable}_previous'].mean()) -
            corr
        )
    else:
        corr, psi_values = np.nan, np.full(len(current_data), np.nan)

    return corr, psi_values

def compute_stats_by_time_age_autocorr(
    df, age_variable, time_indices, age_groups, mean_vars, sd_vars,
    corr_vars, autocorr_vars, regress_vars, time_var='time', member_id_var='member_id'
):
    """
    Compute statistics by time and age group in a DataFrame.

    Calculates mean, standard deviation, pairwise correlations,
    and autocorrelation for specified variables, grouped by time and age.
    Supports half-open or half-closed age groups.

    Parameters:
    ----------
    df : pandas.DataFrame
        DataFrame containing the data.
    age_variable : str
        Column name for age.
    time_indices : list
        Time indices for which statistics are computed.
    age_groups : dict
        Age groups as a dictionary with tuples specifying range.
    mean_vars : list
        Variables for mean calculation.
    sd_vars : list
        Variables for standard deviation calculation.
    corr_vars : list
        Variable pairs for pairwise correlation.
    autocorr_vars : list
        Variables for autocorrelation calculation.
    regress_var : list
        Variables for regression coefficients calculations
    time_var : str, optional
        Column name for time. Default is 'time'.
    member_id_var : str, optional
        Column name for member IDs. Default is 'member_id'.

    Returns
    -------
    tuple
        A tuple containing three elements:
        - DataFrame with computed statistics.
        - Latex regression tables.
        - Dictionary with influence functions for each moment.
        - Dictionary with raw moments.

    Notes
    -----
    - The input DataFrame should be specific to a gender and treatment combination.

    - Output DataFrame (`results_df`):
      This DataFrame uses a MultiIndex for indexing. The first level of the index
      represents age groups, and the second level represents statistical moments 
      (like mean, standard deviation, etc.) calculated for each variable.
      For example, to access the mean of 'var1' for the age group '20-30' at time 0,
      you would use: results_df.loc['20-30', 'mean_var1_t_0'].

    - Output Influence Function Dictionary (`Psi`):
      This dictionary uses a composite key that combines the moment and age group.
      The key format reflects both the statistical moment and the age group it pertains to.
      For instance, the influence function for the mean of 'var1' for the age group '20-30'
      at time 0 can be accessed with the key: Psi['mean_var1_t_0_a_20-30'].
      Here, 'mean_var1_t_0' represents the moment (mean of 'var1' at time 0),
      and 'a_20-30' specifies the age group ('20-30').

    - Output Raw Moments Dictionary:
      Similar to the `Psi` dictionary, this dictionary also uses a composite key format.
      The keys represent both the statistical moment and the age group.
      The structure and access method are identical to that of the `Psi` dictionary.
    """
    results_df = pd.DataFrame()
    Psi = {}
    raw_results = {}
    

    for age_group, (include_min, min_age, max_age, include_max) in age_groups.items():
        age_group_results = {}
        for time in time_indices:
            # Filter DataFrame based on age group and time
            if include_max and include_min:
                age_group_df = df[(df[time_var] == time) & 
                                  (df[age_variable] >= min_age) & 
                                  (df[age_variable] <= max_age)]
            elif include_min:
                age_group_df = df[(df[time_var] == time) & 
                                  (df[age_variable] >= min_age) & 
                                  (df[age_variable] < max_age)]
            else:
                age_group_df = df[(df[time_var] == time) & 
                                  (df[age_variable] > min_age) & 
                                  (df[age_variable] < max_age)]

            # Compute statistics for each variable
            if not age_group_df.empty:
                for var in mean_vars:
                    mean_key = f'mean_{var}_t_{time}_a_{age_group}'
                    age_group_results[f'mean_{var}_t_{time}'] = age_group_df[var].mean()
                    raw_results[mean_key] = age_group_df[var].mean()
                    Psi[mean_key] = age_group_df[var] - age_group_df[var].mean()

                for var in sd_vars:
                    sd_key = f'sd_{var}_t_{time}_a_{age_group}'
                    age_group_results[f'sd_{var}_t_{time}'] = age_group_df[var].std()
                    raw_results[sd_key] = age_group_df[var].std()
                    Psi[sd_key] = (age_group_df[var] - age_group_df[var].mean())**2

                regression_table = pd.DataFrame(columns=['Variable', 'Coefficient'])
                for pair in regress_vars:
                    for var1, var2 in pair.items():
    
                        reg_key = f'reg_{var1}_{var2.replace(", ", "_")}_t_{time}_a_{age_group}'
                        y = age_group_df[var1].values
                        x_vars = var2.split(', ') 
                        x_matrix = age_group_df[x_vars].values
                        X = np.column_stack((np.ones(x_matrix.shape[0]), x_matrix))
                        
                        coefficients = np.linalg.inv(X.T @ X) @ X.T @ y
                        intercept = coefficients[0]
                        slopes = coefficients[1:]

                        # Update the regression table
                        regression_table = pd.concat([
                            regression_table,
                            pd.DataFrame({'Variable': ['Intercept'] + x_vars, 'Coefficient': [intercept] + list(slopes)})
                        ], ignore_index=True)

                        latex_regression_table = regression_table.to_latex(index=False, escape=False)
 
                        age_group_results[f'reg_{var1}_{var2.replace(", ", "_")}_t_{time}'] = str(coefficients)
                        raw_results[reg_key] = str(coefficients)
                        Psi[reg_key] = 0

                for var1, var2 in corr_vars:
                    corr_key = f'corr_{var1}_{var2}_t_{time}_a_{age_group}'
                    age_group_results[f'corr_{var1}_{var2}_t_{time}'] = age_group_df[[var1, var2]].corr().iloc[0, 1]
                    raw_results[corr_key] = age_group_df[[var1, var2]].corr().iloc[0, 1]
                    Psi[corr_key] = ((age_group_df[var1] - age_group_df[var1].mean()) * 
                                     (age_group_df[var2] - age_group_df[var2].mean()))

                # Compute autocorrelation for specified variables
                if time > min(time_indices):
                    for var in autocorr_vars:
                        autocorr_key = f'autocorr_{var}_t_{time}_a_{age_group}'
                        # Custom logic for autocorrelation calculation\
                        current_time = time 
                        previous_time = time_indices[time_indices.index(time) - 1]
                        autocorr, psi = _computeAC(df, var,current_time, previous_time, time_var, member_id_var)
                        age_group_results[f'autocorr_{var}_t_{time}'] = autocorr
                        raw_results[autocorr_key] = autocorr
                        Psi[autocorr_key] = psi

        results_df = pd.concat([results_df, pd.DataFrame(age_group_results, index=[age_group])])

    return results_df, latex_regression_table, Psi, raw_results

def generate_moments(df, config):
    """
    Generate moments for SMM estimation from a DataFrame and a configuration file.
    """

    # Process derived variables under 'Setup'
    derived_vars_config = config['Setup']['derived_variables'][0]

    if derived_vars_config != 'None':
        for derived_var, expression in derived_vars_config.items():
            df[derived_var] = df.eval(expression)

    time_var = config['Setup']['time_var'][0]
    member_id_var = config['Setup']['ID_var'][0]
    treatment_var = config['Setup']['treat_var'][0]

    if treatment_var == 'None':
        treatment_var = False

    gender_var = config['Setup']['gender_var'][0]
    age_var = config['Setup']['age_var'][0]

    # Extract additional setup configurations
    time_indices = config['Setup']['time_indices'][0]

    # Extract configurations under 'Identification'
    mean_vars = config['Identification']['mean']
    sd_vars = config['Identification']['sds']
    corr_vars = [tuple(pair) for pair in config['Identification']['corrs']]
    autocorr_vars = config['Identification']['autocorrs']
    age_groups = {key: tuple(value) for key, value in 
                config['Setup']['age_groups'].items()}
    regress_vars = config['Identification']['OLS']

    # Initialize dictionary for storing DataFrames
    stats_by_group = {}
    raw_moments = {}
    psi_by_group = {}

    # Iterate over combinations of gender and treatment groups
    for gender in df[gender_var].unique():
        
        # Check if treatment_var is specified in the config
        if treatment_var:
            for treatment in df[treatment_var].unique():
                
                df_filtered = df[(df[gender_var] == gender) & 
                                (df[treatment_var] == treatment)]
                stats_df, regression_table, psi_df = compute_stats_by_time_age_autocorr(
                    df_filtered, age_var, time_indices, age_groups, mean_vars, sd_vars,
                    corr_vars, autocorr_vars, regress_vars, time_var=time_var, 
                    member_id_var=member_id_var
                )
                stats_by_group[(gender, treatment)] = stats_df
                psi_by_group = psi_df
        else:
            # Process without filtering by treatment if treatment_var is not specified
            df_filtered = df[df[gender_var] == gender]
            
            stats_df, regression_table, psi_df, age_group_results = compute_stats_by_time_age_autocorr(
                df_filtered, age_var, time_indices, age_groups, mean_vars, sd_vars, 
                corr_vars, autocorr_vars, regress_vars, time_var=time_var, 
                member_id_var=member_id_var
            )
            stats_by_group[(gender, 'None')] = stats_df
            psi_by_group[(gender, 'None')] = psi_df
            raw_moments[(gender, 'None')] = age_group_results

    return stats_by_group, regression_table, psi_by_group, raw_moments


if __name__ == "__main__":
    import yaml
    import pandas as pd

    def read_config(file_path):
        """
        Read a YAML configuration file.
        """
        with open(file_path, 'r') as file:
            return yaml.safe_load(file)
    
    # Reading configuration and data files
    config = read_config('/moments_LS.yml')
    df_in = pd.read_csv('/example_LS.csv')

    # Data format assumptions:
    # - Unique member_ID per individual (row_var)
    # - Columns: member_ID, variables, time, age, gender, treatment (optional)

    # Generate moments and influence functions
    # Output dictionaries for raw moments and influence functions
    # are indexed by gender x treatment group pairs
    # Example: raw_moments[('Female', 'None')] for Female, no treatment group
    df, reg, psi_df, raw_moments = generate_moments(df_in, config)

    # Compute the variance-covariance matrix for the moments
    #V = _compute_cov_matrix_with_nan(raw_moments[('Female', 'None')],
    #                                 psi_df[('Female', 'None')])
