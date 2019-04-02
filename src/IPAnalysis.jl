module IPAnalysis
using Reexport

include("plotting.jl")
@reexport using IPAnalysis.Plotting

include("table_export.jl")
@reexport using IPAnalysis.latex_table

end
