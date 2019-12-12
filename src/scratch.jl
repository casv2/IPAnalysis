using PoSH, JuLIP
using Weave
using Plots

function read_IP(json_file)
    D = load_json(json_file)
    IP = decode_dict(D["IP"])

    Zs = IP.components[2].Vout.zlist.list
    IP_2B = IP.components[2]
    IP_PoSH = IP.components[3]

    return collect(Zs), D["info"], IP_2B, IP_PoSH
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

    plot!(p1, R, energy_2b_11, label="$(s1)-$(s1)")
    plot!(p1, R, energy_2b_12, label="$(s1)-$(s2)")
    plot!(p1, R, energy_2b_22, label="$(s2)-$(s2)")

    @show vcat(energy_2b_11, energy_2b_12, energy_2b_22)

    ymin = minimum(vcat(energy_2b_11, energy_2b_12, energy_2b_22))

    @show ymin

    ylims!(ymin*1.5, -ymin*1.5)
    xlabel!("Interatomic distance [Å]")
    ylabel!("Energy [eV]")

    display(p1)

end

Zs, info, IP_2B, IP_PoSH = read_IP(json_file)

plot_2B(Zs, IP_2B)

#energy_2b(IP_2B, [0,0,2.0], 13, 22)

#function plot_2B(Zs, IP_2B)
