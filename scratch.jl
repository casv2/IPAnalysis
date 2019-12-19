using IPAnalysis

#Zs, info, IP_2B, IP_PoSH, IP = read_IP("TiAl_MD_bulk_liq_CASTEP_all_LD1_2reg_1.361_rep.json")


IP_pdf("TiAl_MD_bulk_liq_CASTEP_all_LD_phonons_1_1_2reg_1.4_rep.json", "2B_PoSH_ph")






# timing(Zs, IP)
#
# plot_2B(Zs, IP_2B)
#
# info["weights"]
#
#
# length(IP_PoSH.coeffs[1])
#
# basis(IP_PoSH)
#
#
#
#
# print("┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n")
# print("┃            WEIGHT TABLE  [E/F/V]         ┃\n")
# print("┣━━━━━━━━━━━━━━━━━┳━━━━━━━━┯━━━━━━━━┯━━━━━━┫\n")
# for cfgs in keys(info["weights"])
#     s = ""
#     if cfgs != "ignore"
#         s *= @sprintf("┃ %15s ┃ %6.2f │ %6.3f │ %4.1f ┃\n", cfgs, info["weights"][cfgs]["E"], info["weights"][cfgs]["F"], info["weights"][cfgs]["V"])
#     end
#     print(s)
# end
# print("┗━━━━━━━━━━━━━━━━━┻━━━━━━━━┷━━━━━━━━┷━━━━━━┛\n")
