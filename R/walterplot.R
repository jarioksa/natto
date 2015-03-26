#' Walter Diagram of Climate
#'
#' Function draws a simplified Walter diagram of monthly temperature
#' and precipitation.
#'
#' Walter diagram displays monthly temperature and precipitation
#' together in the same graph with blue line for precipitation and red
#' line for temperature. The diagram is scaled so that 10 degrees C
#' corresponds to 20 mm precipitation. When precipitation line is
#' above temperature, the period is regarded as humid and coloured
#' skyblue in the diagram, and when precipitation line is below
#' temperature, period is regarded as arid and coloured yellow. Dark
#' blue us used when monthly precipitation exceeds 100 mm.
#'
#' The function draws a simplified diagram intended to be used in
#' general presentations and slides. For proper diagrams you need to
#' find something else.
#'
#' @references Walter, H. (1985) Vegetation of the Earth and
#' ecological systems of the geo-biosphere. 3rd ed. Springer.
#' 
#' @author Jari Oksanen
#'
#' @examples
#' walterplot(c(-9.6,-9.3,-4.8,1.4,7.8,13.5,16.5,14.1,8.9,3.3,-2.8,-7.1),
#' c(31,26,26,20,37,46,71,65,44,45,36,30), "Oulu")
#' 
#' @param temp Vector of twelve monthly temperature values.
#' @param prec Vector of twelve monthly precipitation values.
#' @param title Main title of the graph.
#' @param south (logical): is the location in southern hemisphere?
#' @param ylim limits of y-axis in temperature units. The graph will
#' not be clipped, but the argument can be used to extend the axes.
#' @param mrange Maximum range of y-axis in temperature units.
#' @param lwd Line width of temperature and precipitation lines.
#'
#' @export
`walterplot` <-
    function(temp, prec, title, south = FALSE, ylim, mrange, lwd=3)
{
    cross <- function(x, y, ind)
    {
        x <- x[ind:(ind+1)]
        y <- y[ind:(ind+1)]
        b <- (y[1]-x[1])/(x[2]-x[1]-y[2]+y[1])
        z <- (x[2]-x[1])*b + x[1]
        list(x=b, y=z)
    }
    if (missing(ylim))
        ylim <- 0
    if(is.data.frame(temp))
        temp <- as.vector(as.matrix(temp))
    if(is.data.frame(prec))
        prec <- as.vector(as.matrix(prec))
    if (south) {
        temp <- temp[c(7:12,1:6)]
        prec <- prec[c(7:12,1:6)]
    }
    omn <- tmon <- c(0, 1:12 - 0.5, 12)
    pmon <- tmon
    tmn <- mean(temp[c(1,12)])
    pmn <- mean(prec[c(1,12)])
    tline <- c(tmn, temp, tmn)
    pline <- c(pmn, prec, pmn)
    if (any(pline > 100) && any(pline < 100)) {
        wetrun <- rle(pline > 100)
        beg <- cumsum(wetrun$lengths)
        x <- numeric(length(beg)-1)
        y <- rep(100, length(beg)-1)
        beg <- beg[-length(beg)]
        for (i in 1:length(beg)) {
            z <- cross(rep(100, length(pline)), pline, beg[i])$x
            step <- diff(pmon[beg[i]:(beg[i]+1)])
            x[i] <- z*step + pmon[beg[i]]
        }
        pline <- c(pline, y)
        pmon <- c(pmon, x)
    }
    if (any(pline > 100))
        pline[pline>100] <- (pline[pline>100] - 100)/5 + 100
    pline <- pline/2
    if (any(tline < 0) && any(tline > 0)) {
        zrun <- rle(tline >= 0)
        beg <- cumsum(zrun$lengths)
        x <- numeric(length(beg) -1)
        y <- rep(0, length(beg)-1)
        beg <- beg[-length(beg)]
        for(i in 1:length(beg)) {
            z <- cross(rep(0, length(tline)), tline, beg[i])$x
            x[i] <-  z + pmon[beg[i]]
        }
        tline <- c(tline, y)
        tmon <- c(tmon, x)
    }
    i <- order(pmon)
    pmon <- pmon[i]
    pline <- pline[i]
    i <- order(tmon)
    tmon <- tmon[i]
    tline <- tline[i]
    rng <- range(tline, pline, 0, ylim)
    if (!missing(mrange))
        rng[2] <- max(rng[2], rng[1]+mrange)
    plot(tmon, tline, ylim=rng, xaxs="i",  type="n", ax=FALSE,
         xlab="", ylab="")
    if (any(prec > 100)) {
        wetrun <- rle(pline >= 50)
        end <- cumsum(wetrun$lengths)
        beg <- c(1, end[-length(end)]+1)
        beg <- beg[wetrun$value]
        end <- end[wetrun$value]
        for (i in 1:length(beg))
            polygon(c(pmon[beg[i]],pmon[beg[i]:end[i]],
                      pmon[end[i]]),
                    c(50,pline[beg[i]:end[i]],50),
                    col="midnightblue", border=NA)
       }
    if (all(prec/2 >= temp)) {
        polygon(c(pmon, rev(tmon)),
                c(pmin(pline,50), pmax(rev(tline),0)), col="skyblue",
                border = NA)
    } else if (all(prec/2 <= temp)) {
        polygon(c(pmon, rev(tmon)),
                c(pmin(pline,50), pmax(rev(tline), 0)), col="orange",
                border = NA)
    } else {
        px0 <- 0
        py0 <- pline[1]
        ty0 <- tline[1]
        trun <- rle(c(pmn, prec, pmn)/2 >= c(tmn, temp, tmn))
        end <- cumsum(trun$lengths)
        for (i in 1:length(trun$lengths)) {
            if (i == length(trun$lengths)) {
                x <- x2 <- 12
                y <- pline[length(pline)]
                y2 <- tline[length(tline)] 
            } else {
                cut <- cross(c(tmn, temp, tmn), c(pmn,prec, pmn)/2, end[i])
                step <- diff(omn[end[i]:(end[i]+1)])
                x2 <- x <- cut$x*step + omn[end[i]]
                y2 <- y <- cut$y
            }
            p <- pmon <= x  & pmon >= px0
            t <- tmon <= x & tmon >= px0
            polygon(c(px0, pmon[p], x,x2, rev(tmon[t]),px0),
                    c(py0, pmin(pline[p],50), y,y2, rev(pmax(tline[t],0)), ty0),
                    col = c("orange","skyblue")[as.numeric(trun$value[i])+1],
                    border=NA)
            px0 <- tx0 <- x 
            py0 <- ty0 <- cut$y
        } 
    }
    lines(pmon, pline, col="blue", lwd=lwd, type="l")
    lines(tmon, tline, col="red", lwd=lwd, type="l")
    axis(1, pos=0, at=1:11, labels=rep("", 11))
    abline(v=c(0,12))
    axis(2, at=seq(10*floor(min(tline/10,0)), 10*ceiling(max(tline/10,0)), by=10), labels=FALSE)
    axis(4, at=seq(0, 10*ceiling(max(pline/10)), by=10), labels=FALSE)
    abline(h=0)
}

