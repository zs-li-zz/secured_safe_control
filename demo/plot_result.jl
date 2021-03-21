using PGFPlotsX
latexengine!(PGFPlotsX.PDFLATEX)
time_scale=51
i=3
time_axis=(0:time_scale-1)*0.1
figure = @pgf TikzPicture(
        Axis(
            PlotInc(
            { color => "black", },
            Table([ time_axis, X[i,1:time_scale] ])
                ),
            PlotInc(
            { color => "black", },
            Table([ time_axis, Xkm_hat[i,1:time_scale] ])
                ),
            # PlotInc(
            # { color => "blue", },
            # Table([ time_axis, Xl[i,1:time_scale] ])
            #     ),
            PlotInc(
            { color => "blue", },
            Table([ time_axis, Xls[i,1:time_scale] ])
                )
                )
                )

print_tex(figure)
# PGFPlotsX.pgfsave("cl_no.tex", figure)
