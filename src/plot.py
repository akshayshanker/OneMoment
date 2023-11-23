import pandas as pd
import matplotlib.pyplot as plt

def plot_variables(generated_df, simulated_df, variable_of_interest, variables_to_plot, year, statistic='mean', save_path=None):
    """
    Generate subplots comparing data between generated and simulated datasets for multiple variables of interest.

    Parameters:
    generated_df : dict
        A dictionary containing DataFrames for generated data, where keys are tuples representing gender or ethnicity and other categorical variables.
    simulated_df : dict
        A dictionary containing DataFrames for simulated data, similar to generated_df.
    variable_of_interest : list
        A list of strings specifying the variables of interest, such as gender or ethnicity.
    variables_to_plot : list
        A list of strings representing the variables to be plotted.
    year : int
        An integer representing the year for which data should be plotted.
    statistic : str, optional (default is 'mean')
        A string specifying the statistic to plot, such as 'mean', 'std', etc.
    save_path : str, optional
        A string specifying the file path to save the subplots. If None, the subplots are displayed but not saved.

    Returns:
    None
        The function generates subplots comparing the specified variables between generated and simulated datasets.
    """

    num_plots = len(variable_of_interest)

    # Set up subplots with multiple columns
    fig, axes = plt.subplots(1, num_plots, figsize=(15 * num_plots, 8), sharey=True)

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

        for variable in variables_to_plot:
            full_variable_name = f'{statistic.lower()}_{variable}_t_{year}'

            axes[idx].plot(age_groups, data1[full_variable_name], label=f'Generated data - {full_variable_name}', linestyle='-', marker='s',
                           color='black', linewidth=2)
            axes[idx].plot(age_groups, data2[full_variable_name], label=f'Simulated data - {full_variable_name}', linestyle='--', marker='o',
                           color='black', linewidth=2)

        axes[idx].set_title(f'{variable_of_interest_single}')
        axes[idx].set_xlabel('Age group')

        if idx == 0:
            axes[idx].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol=len(variables_to_plot) * 2)

    axes[-1].set_xlabel('Age group')
    axes[-1].set_ylabel('Variable Values')

    # Save or show subplots
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    else:
        plt.show()