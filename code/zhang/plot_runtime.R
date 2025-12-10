library(ggplot2)
library(patchwork)
options(scipen=999)
source("code/utilities.R")


runtimes <- snakemake@input[['runtimes']]
runtime_plot <- snakemake@output[['plot']]

# read runtimes
cells <- sub(".*_", "", runtimes)
cells <- as.numeric(sub(".txt", "", cells)) * 2 # 2x because we combine fetal and adult cells 1:1
method <- sub("benchmarks/atlas/", "", runtimes)
method <- sub("_.*", "", method)

all.runtimes <- lapply(runtimes, read.table, sep = "\t", header = TRUE)
for (i in seq_along(all.runtimes)) {
    all.runtimes[[i]]$method <- method[i]
    all.runtimes[[i]]$cells <- cells[i]
}

runtimes <- do.call(what = rbind, args = all.runtimes)
runtimes$method <- factor(runtimes$method, levels = c('remo', 'peaks'))

runtimes <- do.call(what = rbind, args = all.runtimes)
runtimes$method <- factor(runtimes$method, levels = c('remo', 'peaks'))
p1 <- ggplot(runtimes, aes(x = cells, y = s/60, color = method, group=method)) + 
  geom_line() +
  geom_point() +
  ylab('Time (minutes)') +
  xlab('Cells') +
  ggtitle("Runtime") +
  scale_color_manual(values = c(REMO_COL, PEAKS_COL)) +
  theme_bw() + theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(200000, 1200000, 200000))+ 
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p2 <- ggplot(runtimes, aes(x = cells, y = max_rss / 1024, color = method, group=method)) + 
  geom_line() +
  geom_point() +
  ylab('Memory (GB)') +
  xlab('Cells') +
  ggtitle("Memory") +
  scale_color_manual(values = c(REMO_COL, PEAKS_COL)) +
  theme_bw() + theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(200000, 1200000, 200000))+ 
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# nonzero elements from each matrix
N_CELL <- seq(50000, 600000, 50000)

matrix_dim <- data.frame()
for (N in N_CELL) {
    fetal_dim <- get_nonzero(paste0("data/fetal/remo_", as.character(N), "/matrix.mtx.gz"))
    adult_dim <- get_nonzero(paste0("data/zhang/remo_", as.character(N), "/matrix.mtx.gz"))
    matrix_dim <- rbind(matrix_dim, data.frame("Nonzero counts" = as.numeric(fetal_dim) + as.numeric(adult_dim),
                                               "Cells" = N*2,
                                               "Method" = "REMO"
                                              ))
    
    fetal_dim <- get_nonzero(paste0("data/fetal/peaks_", as.character(N), "/matrix.mtx.gz"))
    adult_dim <- get_nonzero(paste0("data/zhang/peaks_", as.character(N), "/matrix.mtx.gz"))
    matrix_dim <- rbind(matrix_dim, data.frame("Nonzero counts" = as.numeric(fetal_dim) + as.numeric(adult_dim),
                                               "Cells" = N*2,
                                               "Method" = "Peaks"
                                              ))
}

matrix_dim$Method <- factor(matrix_dim$Method, levels = c("Peaks", "REMO"))

p3 <- ggplot(matrix_dim, aes(x = Cells, y = Nonzero.counts/1e9, color = Method, group=Method)) + 
  geom_line() +
  geom_point() +
  ylab('Nonzero elements (x10^9)') +
  xlab('Cells') +
  ggtitle("Matrix size") +
  geom_hline(yintercept = (2^32)/2/1e9, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c(PEAKS_COL, REMO_COL)) +
  theme_bw() +
  scale_x_continuous(breaks = seq(200000, 1200000, 200000))+ 
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

pp <- (p1 | p2 | p3) + plot_layout(guides = "collect")

ggsave(filename = runtime_plot, plot = pp, height = 3, width = 10, dpi = 500)