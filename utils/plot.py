import pandas as pd
import matplotlib.pyplot as plt

def plot_variables(generated_df, simulated_df, variable_of_interest, variables_to_plot, year, statistic='mean', save_path=None):
    """
    Plots variables from generated and simulated dataframes for a specified variable of interest.

    Parameters:
    - generated_df (pd.DataFrame): DataFrame containing generated data.
    - simulated_df (pd.DataFrame): DataFrame containing simulated data for comparison.
    - variable_of_interest (list): List of variables to be analyzed and plotted.
    - variables_to_plot (list): List of specific variables to be plotted for each age group.
    - year (int): The year for which the analysis is conducted.
    - statistic (str, optional): The type of statistic to use (e.g., 'mean', 'median'). Default is 'mean'.
    - save_path (str, optional): If specified, the plot will be saved to the specified file path. Default is None.

    Returns:
    - None: If save_path is not specified, the plot is displayed. Otherwise, the plot is saved to the specified file path.

    Example:
    >>> gender = ['Female', 'Male']
    >>> variables_to_plot = ['wealth_fin', 'wealth_real']
    >>> year = 2010
    >>> path = 'comparison_plot.png'
    >>> plot_variables(df, modified_data, variable_of_interest = gender, variables_to_plot = variables_to_plot, year = year, save_path = path)
    """
    num_plots = len(variable_of_interest)

    # Set up subplots with multiple columns
    fig, axes = plt.subplots(1, num_plots, figsize=(12 * num_plots, 8), sharey=True)

    for idx, variable_of_interest_single in enumerate(variable_of_interest):
        # Find the first key with the specified variable_of_interest
        found_key_gen = None
        for key_gen in generated_df.keys():
            if variable_of_interest_single.lower() in [v.lower() for v in key_gen]:
                found_key_gen = key_gen
                break

        found_key_sim = None
        for key_sim in simulated_df.keys():
            if variable_of_interest_single.lower() in [v.lower() for v in key_sim]:
                found_key_sim = key_sim
                break

        if not found_key_gen or not found_key_sim:
            print(f"No data found for variable_of_interest: {variable_of_interest_single}")
            continue

        # Extract data for the specified variable_of_interest and variables
        data1 = generated_df[found_key_gen]
        data2 = simulated_df[found_key_sim]

        age_groups = data1.index.tolist()

        # Customize plot
        for variable in variables_to_plot:
            full_variable_name = f'{statistic.lower()}_{variable}_t_{year}'

            axes[idx].plot(age_groups, data1[full_variable_name], label=f'Generated - {variable}', linestyle='-', marker='s',
                           color='black', linewidth=2.5)
            axes[idx].plot(age_groups, data2[full_variable_name], label=f'Simulated - {variable}', linestyle='--', marker='o',
                           color='black', linewidth=2.5)

        # Customize subplot further
        axes[idx].set_title(f'{variable_of_interest_single}')
        axes[idx].set_xlabel('Age cohort')

        # Remove top and right spines
        axes[idx].spines['top'].set_visible(False)
        axes[idx].spines['right'].set_visible(False)

        # Add legend below the first graph
        
        if idx == 0:
            axes[idx].legend(loc='lower center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol=len(variables_to_plot) * 2, handlelength=3)

    # Add labels below the subplots
    axes[-1].set_xlabel('Age cohort')
    # Save or show subplots
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    else:
        plt.show()