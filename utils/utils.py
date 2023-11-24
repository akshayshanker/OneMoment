def create_latex_table(simulated_dict, real_dict, precision=3):
    """
    Creates a LaTeX table comparing simulated and real values with coefficients and standard errors.

    Parameters:
    - simulated_dict (dict): Dictionary containing simulated values for each variable.
    - real_dict (dict): Dictionary containing real values for each variable.
    - precision (int): Number of decimal places to round the coefficients and standard errors (default is 3).

    Returns:
    - latex_table (str): A LaTeX-formatted table as a string, ready to be included in a LaTeX document.
    """
    # Begin LaTeX table
    latex_table = "\\begin{table}[H]\n"
    latex_table += "\\centering \n"
    latex_table += "  \\caption{} \n"
    latex_table += "  \\label{}\n\n"
    latex_table += "\\smallskip\n"
    
    # Begin tabular environment
    num_variables = len(simulated_dict[next(iter(simulated_dict))]['Variable'])
    latex_table += "\\begin{tabular}{@{} l " + " ".join([f"*{{3}}{{d{{5.{precision}}}}}" for _ in range(num_variables * 2)]) + " @{}} \n"
    latex_table += "\\toprule\n"
    
    # Table headers
    variable_names = simulated_dict[next(iter(simulated_dict))]['Variable']
    latex_table += "& \\multicolumn{" + str(num_variables) + "}{c@{}}{Simulated} & \\multicolumn{" + str(num_variables) + "}{c@{}}{Real} \\\\ \n"
    latex_table += "\\cmidrule(lr){" + "2-" + str(num_variables +1) + "}\\cmidrule(lr){"  + str(num_variables +2) + "-" + str(2* num_variables +1) + "}\n"
    latex_table += "& " + " & ".join([f"\\mc{{{name.replace('_', ' ')}}}" for name in variable_names * 2]) + "\\\\ \n"
    latex_table += "\\midrule\n"
    
    # Table content
    for var, values1 in simulated_dict.items():
        values2 = real_dict[var]
        latex_table += f" ${var.replace('_', ' ')}$"
        for i in range(num_variables):
            latex_table += " & " + " & ".join([f"{round(values1['Coefficient'][i], precision)}", f"{round(values2['Coefficient'][i], precision)}"])
        latex_table += "\\\ \n"
        for i in range(num_variables):
            latex_table += " & " + " & ".join([f"({round(values1['Standard Error'][i], precision)})", f"({round(values2['Standard Error'][i], precision)})"])

        
        latex_table += "\\\\ [1ex]\n"

    # End table
    latex_table += "\\bottomrule\n"
    latex_table += "\\end{tabular} \n"
    latex_table += "\\end{table} \n"

    return latex_table