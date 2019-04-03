module latex_table

using JSON
using Printf, Formatting
export to_latex

function to_latex(dict, filename; float_format = "%.3e", cols = [], rows = [], col_names = [], col_code = [])
    if cols == []
        keys_dict = collect(keys(dict))
        nb_col = length(keys(dict))
        col_names = collect(keys(dict))
    else
        keys_dict = cols
        nb_col = length(cols)
        if col_names == []
            col_names = keys_dict
        else
            @assert (length(col_names) == length(keys_dict))
        end
    end
    if rows == []
        rows = collect(1:length(dict[keys_dict[1]]))
    end
    nb_rows = length(rows)
    #Table definition
    latex_table = "\\begin{tabular}{"*repeat("|l", nb_col)*"|} \\hline \n"

    if col_code == []
        for (k,key) in enumerate(keys_dict[1:end-1])
            latex_table *= "\\begin{tabular}{c}" * col_names[k] * "\\end{tabular} & "
        end
        latex_table *= "\\begin{tabular}{c}" * col_names[end] * "\\end{tabular} "
        latex_table *= "\\\\ \\hline \n"
    else
        latex_table *= col_code
    end
    # Write each of the lines
    for i in rows
        for (j,key) in enumerate(keys_dict)
            if typeof(dict[key][i]) == Float64
                # s = sprintf1("%.$(digits)e", dict[key][i])
                s = sprintf1(float_format, dict[key][i])
                latex_table *= s
            else
                latex_table *= string(dict[key][i])
            end
            if j != nb_col
                latex_table *= " & "
            end
        end
        latex_table *=  "\\\\ \\hline \n"
    end
    latex_table *= "\\end{tabular}"

    template = "\\documentclass[a4paper,landscape]{article}
    \\usepackage{booktabs}
    \\usepackage[a4paper,margin=1in,landscape,twocolumn]{geometry}
    \\usepackage{amsmath}
    \\usepackage{graphicx}
    \\usepackage{subfig}
    \\usepackage{diagbox}
    \\usepackage{supertabular}
    \\usepackage{multirow}
    \\usepackage{makecell}
    \\usepackage{hhline}
    \\begin{document}
    \\begin{center}
    \\vspace{3mm}
    $latex_table \\\\
    \\vspace{3mm}
    \\end{center}
    \\end{document}"


    write("out.tex", template)

    filename2 = filename * "_latex_table.pdf"

    run(`pdflatex out.tex`)
    run(`mv out.pdf $filename2`)
    # sleep(1)
    # run(`rm out.tex out.log out.aux`)

end

end
