using NBodyIPs, IPFitting, StaticArrays, Plots, JuLIP, IPAnalysis

IP, info = load_ip("./cvdo_W_4b.json")
IPenv, lsqenf = load_ip("./W_env_r.json")
#IPAnalysis.Plotting.IP_plot(IP)

using NBodyIPs
using Printf, Plots, NBodyIPs, StaticArrays, JuLIP, IPFitting

# function IP_plot(IP::NBodyIPs.NBodyIP; ylim = [-0.8,0.8], xlim = [1.5,8], r0 = 0, return_plot = true, save_plot = false, N = "IP", title = "", filename = "plot.png", θ0 = 0.2)
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
#
# IP_plot(IP, r0=rnn(:W))
#
# IP.components
#
#
# rr = range(2, 8, length=200)
# θ0 = 0.3#acos(-1/3)
#
#
# V2(r) = IP.components[2](r)
# plot(rr, V2.(rr), label="V2")
#
# V3(r1,r2,θ) = IP.components[3]( (SVector(r1, r2), SVector(cos(θ)*r1*r2) ) )
# plot(rr, V3.(rr, rr, θ0), label = "V3 (BA)")
# #plot(rr, V3.(rr,rr,rr), label="V3 (BL)")
#
# V4(r1,r2,r3, θ1, θ2, θ3) = IP.components[4]( (SVector(r1, r2, r3), SVector(θ1, θ2, θ3)) )
# plot(rr, V4.(rr, rr, rr, θ0, θ0, θ0), label="V4 (BA)")
#
# Dict(IP.components[3])["D"]["__id__"]
#
# Dict(IP.components[3])["N"]
#
# Dict(IP.components[2])["Vout"]

function IP_plot(IP, r0; ylims = [-0.1,0.1], xlims = [1.5,8])

    bodyorders = Dict()

    rr = range(xlims[1], xlims[2], length=200)
    θ0 = 0.5

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

    p = plot( yaxis=([ylims[1],ylims[2]]) )

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

    if length(indices_3b) == 1
        try
            V3(r1,r2,r3) = IP.components[indices_3b[1]]( SVector(r1,r2,r3) )
            plot!(p, rr, V3.(rr,rr,rr), label="V3 (BL)")
        catch
            V3(r1,r2,θ) = IP.components[indices_3b[1]]( (SVector(r1, r2), SVector(cos(θ)) ) )
            plot!(p, rr, V3.(rr, rr, θ0), label = "V3 (BA)")
        end
    elseif length(indices_3b) > 1
        for index in indices_3b
            try
                V3(r1,r2,r3) = IP.components[index].Vr( SVector(r1,r2,r3) )
                plot!(p, rr, V3.(rr,rr,rr), label="V3 (BL) (env)")
            catch
                V3(r1,r2,θ) = IP.components[index].Vr( (SVector(r1, r2), SVector(cos(θ)) ) )
                plot!(p, rr, V3.(rr, rr, θ0), label = "V3 (BA) (env)")
            end
        end
    end

    indices_4b = [k for (k,v) in bodyorders if v==4]

    if length(indices_4b) == 1
        try
            V4(r1,r2,r3,r4,r5,r6) = IP.components[indices_4b[1]]( SVector(r1,r2,r3,r4,r5,r6) )
            plot!(p, rr, V4.(rr,rr,rr,rr,rr,rr), label="V4 (BL)")
        catch
            V4(r1,r2,r3, θ1, θ2, θ3) = IP.components[indices_4b[1]]( (SVector(r1, r2, r3), SVector(θ1, θ2, θ3)) )
            plot!(p, rr, V4.(rr, rr, rr, θ0, θ0, θ0), label="V4 (BA)")
        end
    elseif length(indices_4b) > 1
        for index in indices_4b
            try
                V4(r1,r2,r3,r4,r5,r6) = IP.components[index].Vr( SVector(r1,r2,r3,r4,r5,r6) )
                plot!(p, rr, V4.(rr,rr,rr,rr,rr,rr), label="V4 (BL) (env)")
            catch
                V4(r1,r2,r3, θ1, θ2, θ3) = IP.components[index].Vr( (SVector(r1, r2, r3), SVector(θ1, θ2, θ3)) )
                plot!(p, rr, V4.(rr, rr, rr, θ0, θ0, θ0), label="V4 (BA) (env)")
            end
        end
    end

    vline!([r0], label="r0", color="black")
    xlabel!("Interatomic Distance (Angstrom)")
    ylabel!("Energy (eV)")
    display(p)


end


IP_plot(IP, rnn(:W), ylims=[-0.1,0.1])
IP_plot(IPenv, rnn(:W), ylims=[-0.1,0.1])

### FUNCTION ROT PLOT




















#V3(2.89, 2.89, )

r0 = rnn(:W)
rr = range(2, 9, length=200)
θr = range(-1, 1, length=200)

V2(r) = IP.components[2](r)
V3(r1,r2,θ) = IP.components[3]( (SVector(r1, r2), SVector(cos(θ)) ) )
#plot(rr, V3.(rr, rr, θ0), label = "V3 (BA)")
# IP.components[3]
# IPenv.components[3].Vr

θ0 = pi/2
p1 = plot()
plot!(p1, rr, V2.(rr), label="V2")
plot!(p1, rr, V3.(rr, rr, θ0), label = "V3 (BA)")
ylims!(p1, -0.1,0.1)
xlabel!(p1, "Interatomic distance (in A)")
ylabel!(p1, "Energy (in eV)")
vline!(p1, [pi], label="r0", color="black")
title!(p1, "theta = pi/2")

θ0 = -1.0
p2 = plot()
plot!(p2, rr, V2.(rr), label="V2")
plot!(p2, rr, V3.(rr, rr, θ0), label = "V3 (BA)")
ylims!(p2, -0.1,0.1)
xlabel!(p2, "Interatomic distance (in A)")
ylabel!(p2, "Energy (in eV)")
vline!(p2, [pi], label="r0", color="black")
title!(p2, "theta = pi")

plot(p1, p2, layout=(2,1))

savefig("theta_comp.png")

p3 = plot()
V3(r1,r2,θ) = IP.components[3]( (SVector(r1, r2), SVector(cos(θ)*r1*r2) ) )
plot!(p3, θr, V3.(r0, r0, θr), label = "V3 (BA)")
xlabel!(p3, "theta (in rads)")
ylabel!(p3, "Energy (in eV)")
#vline!([pi], label="pi", color="black")

savefig("e_vs_theta.png")


SVector(rr,rr)

IP.components[3](3.0, 3.0, 3.0)

SVector(rr, rr)

SVector(cos(θ0)*rr*rr)

Atoms()
