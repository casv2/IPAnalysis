module Plotting

#using Reexport

#import JuLIP, NBodyIPFitting
using Printf, Plots, NBodyIPs, StaticArrays, JuLIP, IPFitting
using LinearAlgebra

export IP_plot, IP_pdf, Evsθ

function find_bodyorders(IP)
    bodyorders = Dict()

    for (index,component) in enumerate(IP.components)
        try
            E0 = Dict(component)["E0"]
        catch
            try
                bodyorders[index] = Dict(component)["N"]
            catch
                try
                    bodyorders[index] = Dict(component)["Vout"]["N"]
                catch
                    bodyorders[index] = Dict(component)["Vr"]["N"]
                end
            end
        end
    end

    return bodyorders
end

function Evsθ(IP, r0)
    θr = range(-1, 1, length=200)

    bodyorders = find_bodyorders(IP)

    indices_3b = [k for (k,v) in bodyorders if v==3]

    p = plot()

    for index in indices_3b
        if Dict(IP.components[index])["D"]["__id__"] == "BondAngleDesc"
            #if haskey(Dict(IP.components[index]), "Vr")
            try
                V3(r1,r2,θ) = IP.components[index].Vr( (SVector(r1, r2), SVector(cos(θ)) ) )
                plot!(p, θr, V3.(r0, r0, θr), label = "V3 (BA) (env)")
            catch
                V3(r1,r2,θ) = IP.components[index]( (SVector(r1, r2), SVector(cos(θ)) ) )
                plot!(p, θr, V3.(r0, r0, θr), label = "V3 (BA)")
            end
        else
            println("No BA plot")
        end
    end

    xlabel!("Theta (in rads)")
    ylabel!("Energy (eV)")
    display(p)
end

function IP_plot(IP, r0; ylims = [-0.1,0.1], xlims = [1.5,8], θ0=0.3)

    if θ0 < -1 || θ0 > 1
        println("choose θ0 in [-1, 1]")
    end

    bodyorders = find_bodyorders(IP)

    p = plot( yaxis=([ylims[1],ylims[2]]) )
    rr = range(xlims[1], xlims[2], length=200)

    indices_2b = [k for (k,v) in bodyorders if v==2]

    if length(indices_2b) == 1
        V2(r) = IP.components[indices_2b[1]](r)
        plot!(p, rr, V2.(rr), label="V2")
    elseif length(indices_2b) == 2
        V2a(r) = IP.components[indices_2b[1]](r)
        V2b(r) = IP.components[indices_2b[2]](r)
        plot!(p, rr, V2b.(rr) + V2a.(rr), label="V2a + V2b")
    end

    indices_3b = [k for (k,v) in bodyorders if v==3]

    for index in indices_3b
        if haskey(Dict(IP.components[index]), "Vr")
            try
                V3(r1,r2,r3) = IP.components[index].Vr( SVector(r1,r2,r3) )
                plot!(p, rr, V3.(rr,rr,rr), label="V3 (BL) (env)")
            catch
                V3(r1,r2,θ) = IP.components[index].Vr( (SVector(r1, r2), SVector(cos(θ)) ) )
                plot!(p, rr, V3.(rr, rr, θ0), label = "V3 (BA) (env)")
            end
        else
            try
                V3(r1,r2,r3) = IP.components[index]( SVector(r1,r2,r3) )
                plot!(p, rr, V3.(rr,rr,rr), label="V3 (BL)")
            catch
                V3(r1,r2,θ) = IP.components[index]( (SVector(r1, r2), SVector(cos(θ)) ) )
                plot!(p, rr, V3.(rr, rr, θ0), label = "V3 (BA)")
            end
        end
    end

    indices_4b = [k for (k,v) in bodyorders if v==4]

    for index in indices_4b
        if haskey(Dict(IP.components[index]), "Vr")
            try
                V4(r1,r2,r3,r4,r5,r6) = IP.components[index].Vr( SVector(r1,r2,r3,r4,r5,r6) )
                plot!(p, rr, V4.(rr,rr,rr,rr,rr,rr), label="V4 (BL) (env)")
            catch
                V4(r1,r2,r3, θ1, θ2, θ3) = IP.components[index].Vr( (SVector(r1, r2, r3), SVector(θ1, θ2, θ3)) )
                plot!(p, rr, V4.(rr, rr, rr, θ0, θ0, θ0), label="V4 (BA) (env)")
            end
        else
            try
                V4(r1,r2,r3,r4,r5,r6) = IP.components[index]( SVector(r1,r2,r3,r4,r5,r6) )
                plot!(p, rr, V4.(rr,rr,rr,rr,rr,rr), label="V4 (BL)")
            catch
                V4(r1,r2,r3, θ1, θ2, θ3) = IP.components[index]( (SVector(r1, r2, r3), SVector(θ1, θ2, θ3)) )
                plot!(p, rr, V4.(rr, rr, rr, θ0, θ0, θ0), label="V4 (BA)")
            end
        end
    end

    vline!([r0], label="r0", color="black")
    xlabel!("Interatomic Distance (Angstrom)")
    ylabel!("Energy (eV)")
    display(p)

end


#export IP_plot



# gr(size=(800,500), html_output_format=:png)
#
# # function unfold(A)
# #     V = []
# #     for x in A
# #         if x === A
# #             push!(V, x)
# #         else
# #             append!(V, unfold(x))
# #         end
# #     end
# #     V
# # end
#
# using NBodyIPs
# using Printf, Plots, NBodyIPs, StaticArrays, JuLIP, IPFitting
#
# function IP_plot(IP::NBodyIPs.NBodyIP; ylim = [-0.8,0.8], xlim = [1.5,8], r0 = 0, return_plot = true, save_plot = false, N = "IP", title = "", filename = "plot.png", θ0 = 0.3)
#     # collect the IPs
#     IPs = NBodyIPs.bodyorder.(IP.components)[2:end]
#
#     rr = range(xlim[1], xlim[2], length=200)
#     #θ0 = 0.3#acos(-1/3)
#
#     p = plot( yaxis=([ylim[1],ylim[2]]) )
#
#     j = 2
#
#     # plot V2a + V2b or V2
#
#     if length(findall(IPs .== 2)) > 1
#         # V2a(r) = IP.components[j](r)
#         j += 1
#         # V2b(r) = IP.components[j](r)
#         # plot!(p, rr, V2b.(rr) + V2a.(rr), label="V2a + V2b")
#         plot!(p, rr, IP.components[2].(rr)+IP.components[3].(rr),  label="V2a + V2b")
#     else
#         V2(r) = IP.components[j](r)
#         plot!(p, rr, V2.(rr), label="V2")
#     end
#
#     j += 1
#
#     # plot 3/4 BA/BL
#
#     for i in deleteat!(IPs, findall((in)([2]), IPs))
#         if i == 3
#             try
#                 V3(r1,r2,r3) = IP.components[j]( SVector(r1,r2,r3) )
#                 plot!(p, rr, V3.(rr,rr,rr), label="V3 (BL)")
#             catch
#                 V3(r1,r2,θ) = IP.components[j]( (SVector(r1, r2), SVector(cos(θ)*r1*r2) ) )
#                 plot!(p, rr, V3.(rr, rr, θ0), label = "V3 (BA)")
#             end
#         elseif i == 4
#             try
#                 V4(r1,r2,r3,r4,r5,r6) = IP.components[j]( SVector(r1,r2,r3,r4,r5,r6) )
#                 plot!(p, rr, V4.(rr,rr,rr,rr,rr,rr), label="V4 (BL)")
#             catch
#                 V4(r1,r2,r3, θ1, θ2, θ3) = IP.components[j]( (SVector(r1, r2, r3), SVector(θ1, θ2, θ3)) )
#                 plot!(p, rr, V4.(rr, rr, rr, θ0, θ0, θ0), label="V4 (BA)")
#             end
#         end
#         j += 1
#     end
#
#     # add r0
#
#     if typeof(r0) == Float64
#         vline!([r0], label="r0", color="black")
#     end
#
#     xlabel!("Interatomic distance (Angstrom)")
#     ylabel!("Energy (eV)")
#     title!(title)
#
#     if save_plot
#         savefig(filename)
#     end
#     if return_plot
#         return p
#     end
# end

function IP_pdf(IP::NBodyIPs.NBodyIP, info::Dict{String,Any}, filename)
    #IP_plot(IP, save_plot = true, filename=filename)
    #error table
    error_table = "\\begin{supertabular}{ l c c c } \\toprule \n"
    error_table *= "Config type & E (eV) & F (eV/A) & V (eV/A2) \\\\ \\midrule \n"

    types = sort(collect(keys(info["errors"]["rmse"])))

    for key in deleteat!(types, findall((in)(types, ["set"])))
        E = try string(info["errors"]["rmse"][key]["E"])[1:7] catch; "NaN" end
        F = try string(info["errors"]["rmse"][key]["F"])[1:7] catch; "NaN" end
        V = try string(info["errors"]["rmse"][key]["V"])[1:7] catch; "NaN" end
        s = @sprintf "%s & %s & %s & %s \\\\ \n" replace(key, "_" => "\\_") E F V
        error_table *= s
    end

    E = try string(info["errors"]["rmse"]["set"]["E"])[1:7] catch; "NaN" end
    F = try string(info["errors"]["rmse"]["set"]["F"])[1:7] catch; "NaN" end
    V = try string(info["errors"]["rmse"]["set"]["V"])[1:7] catch; "NaN" end
    s = @sprintf "%s & %s & %s & %s \\\\ \n" "set" E F V
    error_table *= s

    error_table *= "\\end{supertabular}"

    #dataweights
    data_table = "\\begin{supertabular}{ l c c c } \\toprule \n"
    s = @sprintf "Data & %s & %s & %s \\\\ \\midrule \n" "E" "F" "V"
    data_table *= s
    s = @sprintf "Weight & %s & %s & %s \\\\ \\midrule \n" info["obsweights"]["E"] info["obsweights"]["F"] info["obsweights"]["V"]
    data_table *= s
    data_table *= "\\end{supertabular}"

    #weighttable
    weight_table = "\\begin{supertabular}{ l c } \\toprule \n"
    weight_table *= "Config type & Weight \\\\ \\midrule \n"

    for key in sort(collect(keys(info["configweights"])))
        s = @sprintf "%s & %s \\\\ \n" replace(key, "_" => "\\_") info["configweights"][key]
        weight_table *= s
    end
    weight_table *= "\\end{supertabular}"

    reg_table = "\\begin{supertabular}{ l c c c} \\toprule \n"
    reg_table *= "Type & r0 & r1 & creg \\\\ \\midrule \n"

    for i in 1:length(info["regularisers"])
        s = @sprintf "%s & %s & %s & %s \\\\ \n" info["regularisers"][i]["type"] info["regularisers"][i]["r0"] info["regularisers"][i]["r1"] info["regularisers"][i]["creg"]
        reg_table *= s
    end

    reg_table *= "\\end{supertabular}"

    db = replace(info["dbpath"][3:end], "_" => "\\_")

    lname = replace(filename, "_" => "\\_")
    #pname = @sprintf("%s.png", filename)

    e0 = info["E0"]

    if e0 == nothing
        e0 = "2B"
    end

    basis = string(length(info["Ibasis"]))

    @show lname, db, e0, basis, data_table, weight_table, error_table, reg_table #, pname

    template = "\\documentclass[a4paper,landscape]{article}
    \\usepackage{booktabs}
    \\usepackage[a4paper,margin=1in,landscape,twocolumn]{geometry}
    \\usepackage{amsmath}
    \\usepackage{graphicx}
    \\usepackage{subfig}
    \\usepackage{diagbox}
    \\usepackage{supertabular}
    \\begin{document}
    \\begin{center}
    \\textbf{Name}: $lname \\\\
    \\vspace{3mm}
    \\textbf{LsqDB}: $db \\\\
    \\vspace{2mm}
    \\textbf{E0}: $e0 \\\\
    \\vspace{2mm}
    \\textbf{Basis functions}: $basis \\\\
    \\vspace{3mm}
    $data_table \\\\
    \\vspace{3mm}
    $reg_table \\\\
    \\vspace{3mm}
    $error_table \\\\
    \\vspace{3mm}
    $weight_table \\\\
    \\end{center}
    \\end{document}"


    write("out.tex", template)

    filename2 = filename * "_IPanalysis.pdf"

    run(`pdflatex out.tex`)
    run(`mv out.pdf $filename2`)
    sleep(1)
    run(`rm out.tex out.log out.aux`)
end

end # module


# \\vspace{3mm}
# \\begin{figure}[h]
#     \\centering
#     \\subfloat{{\\includegraphics[height=8cm]{$pname} }}%
#     \\caption{Slices of \$V_{n}\$}%
# \\end{figure}

# function force_plot(IP::NBodyIPs.NBodyIP, test_data::Array{NBodyIPFitting.Dat,1}; s = 10, return_plot = true)
#     data = sortcols(hcat([[split(test_data[i].configtype, ":")[1], test_data[i].D["F"], unfold(forces(IP, JuLIP.Atoms(test_data[i])))] for i in  1:s:length(test_data)]...))
#     uniq_config_sel = unique(data[1,:])
#
#     p1 = plot()
#     p2 = plot()
#
#     ylabel!(p1, "Predicted Force (eV/A)")
#     xlabel!(p1, "Target Force (eV/A)")
#     title!(p1, "Force Scatter Plot")
#
#     ylabel!(p2, "Predicted Force Error (eV/A)")
#     xlabel!(p2, "Target Force (eV/A)")
#     title!(p2, "Force Error Plot")
#
#     plot_max = findmax(unfold(data[3,:]))[1]
#     plot_min = findmin(unfold(data[3,:]))[1]
#     for uniq_config in uniq_config_sel
#         indices = find(data[1,:] .== uniq_config)
#         min, max = findmin(indices)[1], findmax(indices)[1]
#         p1 = scatter!(p1, unfold(data[2, min:max]), unfold(data[3, min:max]), label=uniq_config, legend=:bottomright)
#         p2 = scatter!(p2, unfold(data[2, min:max]), abs.(unfold(data[3, min:max]) - unfold(data[2, min:max])), label=uniq_config, legend=false, yaxis=(:log10, (0.00005,2)))
#     end
#
#     P = plot!(p1, [plot_min, plot_max], [plot_min, plot_max], color="black")
#
#     P = plot(p1, p2, layout=(1,2), size=(1400,700))
#
#     if return_plot == false
#         savefig("forceplot.png")
#     else
#         display(P)
#     end
# end
#
# function energy_plot(IP::NBodyIPs.NBodyIP, test_data::Array{NBodyIPFitting.Dat,1}; s = 10, return_plot = true)
#     data = sortcols(hcat([[split(test_data[i].configtype, ":")[1], (test_data[i].D["E"][1]/length(test_data[i].at.Z)), (energy(IP, JuLIP.Atoms(test_data[i]))/length(test_data[i].at.Z))] for i in  1:s:length(test_data)]...))
#
#     uniq_config_sel = unique(data[1,:])
#
#     p1 = plot()
#     p2 = plot()
#
#     ylabel!(p1, "Predicted Energy per Atom (eV)")
#     xlabel!(p1, "Target Energy per Atom(eV)")
#     title!(p1, "Energy Scatter Plot")
#
#     ylabel!(p2, "Predicted Energy per Atom (eV)")
#     xlabel!(p2, "Target Energy per Atom (eV)")
#     title!(p2, "Energy Error Plot")
#
#     plot_max = findmax(unfold(data[3,:]))[1]
#     plot_min = findmin(unfold(data[3,:]))[1]
#
#     for uniq_config in uniq_config_sel
#         indices = find(data[1,:] .== uniq_config)
#         min, max = findmin(indices)[1], findmax(indices)[1]
#         p1 = scatter!(p1, data[2, min:max], data[3, min:max], label=uniq_config, legend=:bottomright)
#         p2 = scatter!(p2, data[2, min:max], abs.(data[3, min:max] - data[2, min:max]), label=uniq_config, legend=false, yaxis=(:log10, (0.00005,0.02)))
#     end
#
#     p1 = plot!(p1, [plot_min, plot_max], [plot_min, plot_max], color="black")
#
#     P = plot(p1, p2, layout=(1,2), size=(1400,700))
#
#     if return_plot == false
#         savefig("energyplot.png")
#     else
#         display(P)
#     end
#
# end

#\\begin{figure}[h]
#    \\centering
#    \\subfloat{{\\includegraphics[height=8cm]{IP_plot.png} }}%
#    \\caption{Slices of \$V_{n}\$}%
#\\end{figure}
