 Pertubation_graph <- plot_ly(data.frame(temp), x = ~Age, y = ~res,
         main = paste( "Pig ID:", i, ", Lambda =", lambda, "\nDifference between CFI and TTC"),
         xlab = "Age (days)",
         ylab = "Amount of difference: CFI - TTC (kg)",
         type = "p", pch = 10, cex = 0.5,
         ylim = c(min(B$dif.CFI, res),
                  max(B$dif.CFI, res)),
         cex.main = 1.5, cex.lab = 1.2, axes = F)
    axis(1, at= seq(dif$eval_day[1],dif$eval_day[length(dif$eval_day)], by = 5), cex.axis = 1.1)
    axis(2, at=, cex.axis = 1.1)
    axis(4, at=, cex.axis = 1.1, col.axis = "blue", col = "blue")
    mtext("Percentage of difference: CFI - TTC (%)", side=4, line = 3, cex = 1.2, col= "blue")
    box()
    abline(crit1,0, col = "red", lty = 2)
    abline(0,0, col = "red")
 saveWidget(ggplotly(Pertubation_graph), file=paste0("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box/Graphs/Step3_graphs/", idc, ".", ID[idc], ".html"))