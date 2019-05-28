#!/usr/bin/env Rscript
library(IRanges)
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(grid)
library(patchwork)
library(ggpubr)
library(thacklr)


# color scale only
max_evalue <- 1
min_evalue <- 1e-15

evalue_color <- function(evalue, min_evalue, max_evalue){
  nle <- -log10(evalue)
  min_nle <- -log10(min_evalue)
  max_nle <- -log10(max_evalue)
  sapply(nle, function(e) max(max_nle, min(e, min_nle)))
}

## setwd("../dev-data/hhplot")
## args <- c("bara_KP-head.tsv", "barHuA8-body.tsv", "foo.pdf", "../../data/")
## grep -P '^(ACC|DESC)' ~/mnt/mit/nobackup1/chisholmlab/databases/pfam/Pfam-32.0/Pfam-A.hmm | sed 's/\s\+/\t/' > pfam-acc-desc.tsv
## pfam_0 <- read_tsv("pfam-acc-desc.tsv", col_names=c("key", "value")) %>%
##   mutate(group=ceiling(row_number()/2)) %>%
##   spread(key, value) %>%
##   select(pfam_acc=ACC, description=DESC)
## write_rds(pfam_0, "hhplot-pfam-32.0.rds")

args <- commandArgs(trailingOnly=TRUE)
head_0 <- read_tsv(args[1], col_names=FALSE) %>% deframe()
hits_0 <- read_tsv(args[2]) %>% mutate(i=row_number())

if(length(args) == 4){
  write("trying to match accessions to pfam descriptions", stderr())
  pfam_1 <- read_rds(paste0(args[4], "/hhplot-pfam-32.0.rds"))
  hits_0 <- hits_0 %>%
    mutate(pfam_acc = str_replace(description, " .*", ""), description = NULL) %>%
    left_join(pfam_1) %>%
    select(i, description, pfam_acc, everything())
}

xmax <- as.numeric(head_0["Match_columns"])
top_rows <- 1:20

hitr <- IRanges(start=hits_0$query_start, end=hits_0$query_end)

bin_size <- 20
binr <- IRanges(start=seq(1, xmax, bin_size), end=seq(1, xmax, bin_size)+bin_size-1)

hits_1 <- findOverlaps(hitr, binr) %>% as_tibble %>%
  rename(i=queryHits, bin=subjectHits) %>%
  right_join(hits_0)

hits_2 <- hits_1 %>% group_by(bin) %>% top_n(10, -evalue) %>%
  ungroup %>% select(-bin) %>% unique %>% arrange(evalue) %>%
  mutate(y=-row_number())

hits_2 %<>% mutate(
  target_full_start = query_start - target_start,
  target_full_end = query_end + (target_cols - target_end)
)


gg_hits <- ggplot(hits_2) +
#  geom_segment(aes(x=(query_start+query_end)/2, xend=xmax*1.08, y=y, yend=(y+1)*2-1), hits_0[top_rows,], alpha=.2) +
  geom_vline(xintercept = c(0,xmax), linetype=3) +
  geom_vline(xintercept = c(-0.05*xmax,1.05*xmax)) +
  geom_segment(aes(x=target_full_start, xend=target_full_end, y=y, yend=y),
               color="black", size=1) +
  geom_segment(aes(x=query_start, xend=query_end, y=y, yend=y,
                   color=evalue_color(evalue, min_evalue, max_evalue)), size=3) + #-log10(evalue))) +
  geom_point(aes(x=query_end*1.04, y=y+.1), color="grey90", size=5) +
  geom_text(aes(x=query_end*1.042, y=y, label=i), size=4, hjust=0.5) +
  #  geom_text(aes(y=(y+1)*2-1, label=description), x=xmax*1.1, hits_0[top_rows,], size=3, hjust=0) +
  expand_limits(color=-log10(c(max_evalue, min_evalue))) +
  scale_color_distiller(palette="Spectral") +
  ggtitle(str_replace(args[3], ".pdf", "")) +
  coord_cartesian(xlim=c(0, xmax), ylim=c(-50,0)) +
  theme_classic() +
  theme(legend.position = "bottom") + no_y_axis()

gg_table <- ggtexttable(hits_2[1:30,1:5], rows = NULL)

gg <- gg_hits + gg_table + plot_layout(widths=c(1,1))
ggsave(args[3], gg, width=16, height=11)
