## Other functions for working with APT ##

plot.tempTraj <- function(apt) {
    traj <- apt$tempTraj
    nRow <- nrow(traj)
    nCol <- ncol(traj)
    myCols <- rainbow(nCol)
    par(mfrow=n2mfrow(2))
    YLIM <- range(pretty(c(0, apt$tempTraj[,nCol])))
    plot(1:nRow, apt$tempTraj[,nCol], typ="n", ylim=YLIM, xlab="Iteration", ylab="Temperature")
    for (ii in 1:nCol)
        lines(apt$tempTraj[,ii], col=myCols[ii])
    YLIM <- range(pretty(c(0, log10(apt$tempTraj[,nCol]))))
    plot(1:nRow, log10(apt$tempTraj[,nCol]), typ="n", ylim=YLIM, xlab="Iteration", ylab="log10(Temperature)")
    for (ii in 1:nCol)
        lines(log10(apt$tempTraj[,ii]), col=myCols[ii])
}
