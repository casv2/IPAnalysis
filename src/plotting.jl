module Plotting

using PoSH, JuLIP
using Plots, Weave, Printf

export read_IP, energy_2b, plot_2B, IP_pdf, timing, basis, weights

function read_IP(json_file)
    D = load_json(json_file)
    IP = decode_dict(D["IP"])

    Zs = IP.components[2].Vout.zlist.list
    IP_2B = IP.components[2]
    IP_PoSH = IP.components[3]

    return collect(Zs), D["info"], IP_2B, IP_PoSH, IP
end

function energy_2b(IP_2B, v, z0, Z)
    Rs = [JVecF(v)]
    z0 = z0
    Zs = Int16[ Z ]
    tmp = JuLIP.alloc_temp(IP_2B, length(Rs))
    return convert(Float64,JuLIP.Potentials.evaluate!(tmp, IP_2B, Rs, Zs, z0))
end

function plot_2B(Zs, IP_2B; xmin = 1, xmax = 6)
    R = collect(range(xmin,xmax,length=100))
    _R = [ [0.0, 0.0, i] for i in range(xmin, xmax, length=100) ]

    p1 = plot()

    if length(Zs) == 2
        z1, z2 = Zs[1], Zs[2]
        s1, s2 = chemical_symbol(z1), chemical_symbol(z2)
        energy_2b_11 = [ energy_2b(IP_2B, r, z1, z1) for r in _R]
        energy_2b_12 = [ energy_2b(IP_2B, r, z1, z2) for r in _R]
        energy_2b_22 = [ energy_2b(IP_2B, r, z2, z2) for r in _R]
    end

    vline!(p1, [rnn(Symbol(s1))], label="rnn $(s1)", color="grey", linestyle=:dash)
    vline!(p1, [rnn(Symbol(s2))], label="rnn $(s2)", color="black", linestyle=:dash)
    plot!(p1, R, energy_2b_11, label="$(s1)-$(s1)")
    plot!(p1, R, energy_2b_12, label="$(s1)-$(s2)")
    plot!(p1, R, energy_2b_22, label="$(s2)-$(s2)")

    ymin = minimum(vcat(energy_2b_11, energy_2b_12, energy_2b_22))

    ylims!(ymin*1.5, -ymin*1.5)
    xlabel!("Interatomic distance [Å]")
    ylabel!("Energy [eV]")

    display(p1)
end

function timing(Zs, IP)
    at = bulk(Symbol(chemical_symbol(Zs[1]))) * (2,2,2)

    for i in 1:3
        at = rattle!(at,0.01)
        forces(IP, at)
    end

    N = 10

    start = time()
    for i in 1:N
        at = rattle!(at,0.01)
        forces(IP, at)
    end
    elapsed = time() - start

    T =  ((elapsed/length(at))/N) * 1000
    T = round(T, digits=3)

    println("TIMING: $(T) [ms/atom]")
end

function basis(IP_PoSH)
    B = length(IP_PoSH.coeffs[1])

    println("#BASIS: $(B)")
end


function weights(info)
    print("┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n")
    print("┃            WEIGHT TABLE  [E/F/V]         ┃\n")
    print("┣━━━━━━━━━━━━━━━━━┳━━━━━━━━┯━━━━━━━━┯━━━━━━┫\n")
    for cfgs in keys(info["weights"])
        s = ""
        if cfgs != "ignore"
            s *= @sprintf("┃ %15s ┃ %6.2f │ %6.3f │ %4.1f ┃\n", cfgs, info["weights"][cfgs]["E"], info["weights"][cfgs]["F"], info["weights"][cfgs]["V"])
        end
        print(s)
    end
    print("┗━━━━━━━━━━━━━━━━━┻━━━━━━━━┷━━━━━━━━┷━━━━━━┛\n")
end

function IP_pdf(json_file, model_name)

    template =
"---
title : $(model_name) report
---

```{julia;echo=false}
using PoSH, IPAnalysis, IPFitting

Zs, info, IP_2B, IP_PoSH, IP = read_IP(\"$(json_file)\")


timing(Zs, IP)

weights(info)

rmse_table(info[\"errors\"])

plot_2B(Zs, IP_2B)
```"

    write("$(model_name).jmd", template)

    weave("$(model_name).jmd", out_path = :pwd, doctype = "pandoc2pdf")

end

end
