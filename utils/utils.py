def export_latex_table(results_dict1, results_dict2):
    print("\\begin{tabular}{lrr|lrr}")
    print("\\toprule")
    print("\\multicolumn{3}{c}{Data} & \\multicolumn{3}{c}{Simulated} \\\\")
    print("\\cmidrule(lr){1-3} \\cmidrule(lr){4-6}")
    print("Variable & Coefficient & Standard Error & Variable & Coefficient & Standard Error \\\\")

    for variable, data1, data2 in zip(results_dict1.keys(), results_dict1.values(), results_dict2.values()):
        formatted_variable = variable.replace('_', ' ')
        print("\\midrule")
        print(f"{formatted_variable} \\\\")
        print("\\hline")
        for var, coef, std_err in zip(data1['Variable'], data1['Coefficient'], data1['Standard Error']):
            # Replace lowercase letters with spaces in the variable name
            formatted_var = var.replace('_', ' ')
            coef_str = f"${coef:.9g}$"  # Format coefficient to scientific notation with 12 significant figures
            std_err_str = f"${std_err:.9g}$"  # Format standard error to scientific notation with 12 significant figures
            print(f"{formatted_var} & {coef_str} & {std_err_str} & {formatted_var} & {data2['Coefficient'][data2['Variable'].index(var)]:.12g} & {data2['Standard Error'][data2['Variable'].index(var)]:.12g} \\\\")

    print("\\bottomrule")
    print("\\end{tabular}")