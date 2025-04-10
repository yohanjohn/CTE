using CairoMakie, Random, LinearAlgebra
import XLSX
fname = "SGI_Raw_Data_Yohan.xlsx"
#data = XLSX.readxlsx("SGI_Raw_Data_Yohan.xlsx")

ESsupra = Float64.(XLSX.readdata(fname,"Early Stage 3","C2:C11"))
ESinfra = Float64.(XLSX.readdata(fname,"Early Stage 3","D2:D11"))
ESSGI = Float64.(XLSX.readdata(fname,"Early Stage 3","F2:F11"))

LSsupra = Float64.(XLSX.readdata(fname,"Late Stage 3","C2:C11"))
LSinfra = Float64.(XLSX.readdata(fname,"Late Stage 3","D2:D11"))
LSSGI = Float64.(XLSX.readdata(fname,"Late Stage 3","F2:F11"))
fig = Figure()

ax1 = Axis(fig[1, 1], title = "A")

#ax2 = Axis(fig[2, 1], title = "B")
#ax3 = Axis(fig[2, 1], title = "C")
#ax4 = Axis(fig[2, 2], title = "D")
# lines!(ax1,1:length(ESsupra),ESsupra[:])
# lines!(ax1,1:length(ESsupra),LSsupra[:])
# lines!(ax2,1:length(ESsupra),ESinfra[:])
# lines!(ax2,1:length(ESsupra),LSinfra[:])
lines!(ax1,1:length(ESSGI),ESSGI[:])
lines!(ax1,1:length(LSSGI),LSSGI[:])
fig
