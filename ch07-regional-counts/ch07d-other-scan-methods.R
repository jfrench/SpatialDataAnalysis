set.seed(112)
library(smerc)
library(sf)

data("neast")
data("neastw")

# setup up arguments
coords <-  neast[, c("x", "y"), drop = TRUE]
pop <- neast$population
cases <- neast$cases

plot(neast["cases"], pal = hcl.colors)

# run various tests

# precog
out1 <- precog.test(coords = coords, cases = cases,
                    pop = pop, w = neastw, nsim = 999)
# circular
out2 <- scan.test(coords = coords, cases = cases,
                  pop = pop, nsim = 999)

# elliptic
out3 <- elliptic.test(coords = coords, cases = cases,
                      pop = pop, nsim = 999, a = 0)

# flex
out4 <- flex_test(coords = coords, cases = cases, pop = pop,
                  w = neastw, k = 15, nsim = 999)

# rflex
out5 <- rflex.test(coords = coords, cases = cases,
                   pop = pop, w= neastw, k = 50, nsim = 999)

# dc
out6 <- dc.test(coords = coords, cases = cases, pop = pop,
                w = neastw, nsim = 999)

# fast
out7 <- fast.test(coords = coords, cases = cases, pop = pop,
                  nsim = 999)

# uls
out8 <- uls.test(coords = coords, cases = cases, pop = pop,
                 w = neastw, nsim = 999)

# create plots of results

# setup plotting arguments
# color palette
col_pal <-  c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
              "#66A61E", "#E6AB02", "#A6761D", "#666666",
              "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
              "#66A61E", "#E6AB02", "#A6761D", "#666666")
# number of detected clusters
nc1 <- nclusters(out1)
# color clusters for out1
col1 <- color.clusters(out1, col = col_pal[seq_len(nc1)])
nc2 <- nclusters(out2)
col2 <- color.clusters(out2, col = col_pal[seq_len(nc2)])
nc3 <- nclusters(out3)
col3 <- color.clusters(out3, col = col_pal[seq_len(nc3)])
nc4 <- nclusters(out4)
col4 <- color.clusters(out4, col = col_pal[seq_len(nc4)])
nc5 <- nclusters(out5)
col5 <- color.clusters(out5, col = col_pal[seq_len(nc5)])
nc6 <- nclusters(out6)
col6 <- color.clusters(out6, col = col_pal[seq_len(nc6)])
nc7 <- nclusters(out7)
col7 <- color.clusters(out7, col = col_pal[seq_len(nc7)])
nc8 <- nclusters(out8)
col8 <- color.clusters(out8, col = col_pal[seq_len(nc8)])

# precog
plot(st_geometry(neast), col = col1)
title("(i) PreCoG", cex.main = 2, line = -1)

# circular
plot(st_geometry(neast), col = col2)
title("(ii) circular", cex.main = 2, line = -1)

# elliptic
plot(st_geometry(neast), col = col3)
title("(iii) elliptic", cex.main = 2, line = -1)

# flex
plot(st_geometry(neast), col = col4)
title("(iv) flex15", cex.main = 2, line = -1)
# rflex
plot(st_geometry(neast), col = col5)
title("(v) rflex", cex.main = 2, line = -1)

# dc
plot(st_geometry(neast), col = col6)
title("(vi) dc", cex.main = 2, line = -1)

# fast
plot(st_geometry(neast), col = col7)
title("(vii) fast", cex.main = 2, line = -1)

# uls
plot(st_geometry(neast), col = col8)
title("(viii) uls", cex.main = 2, line = -1)


# table comparing methods

# function to get relevant data
get_table_data <- function(x, method) {
  data.frame(
    method = method,
    cluster = paste("cluster", seq_along(x$clusters)),
    pop = sgetElement(x$clusters, "population"),
    cases = sgetElement(x$clusters, "cases"),
    ex = sgetElement(x$clusters, "expected"),
    smr = sgetElement(x$clusters, "smr"),
    nregions = sapply(lgetElement(x$clusters, "locids"), length),
    pvalue = sgetElement(x$clusters, "pvalue")
  )
}

# data in list format
table_list <- mapply(get_table_data,
                     x = list(out1, out2, out3, out4,
                              out5, out6, out7, out8),
                     method = c("precog", "circular",
                                "elliptic", "flex15",
                                "rflex", "dc", "fast",
                                "uls"),
                     SIMPLIFY = FALSE)
# create data frame
table_data <- do.call(rbind, table_list)
# format output
table_data$ex <- round(table_data$ex)
table_data$smr <- round(table_data$smr, 3)
# final results
table_data
