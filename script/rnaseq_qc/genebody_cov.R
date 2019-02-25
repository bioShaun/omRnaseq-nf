options(warn = -1)
options(stringsAsFactors = F)
suppressMessages(library('reshape2'))
suppressMessages(library('argparser'))
suppressMessages(library('plyr'))
suppressMessages(library('dplyr'))
suppressMessages(library('omplotr'))

p <- arg_parser('genebody cov')
p <- add_argument(p,'--cov_dir',help = 'genebody cov directory.')
argv <- parse_args(parser = p)


gene_cov_plot <- function(filter_df, output) {
  sample_number <- length(unique(filter_df$sample))
  p <- ggplot(filter_df, aes(percentile, reads_portion, color=sample)) +
    geom_point() +
    geom_line() +
    theme_onmath() +
    scale_y_continuous(labels = scales::percent) +
    xlab('Percentile') + ylab('Reads Portion') +
    scale_color_brewer(palette = 'Set1')
  
  if (sample_number > 1) {
    facet_wrap_ncol = round(sqrt(sample_number))
    p <- p + guides(color=guide_legend(title = ''))
  } else {
    p <- p + guides(color=F)
  }
  
  if (! is.null(output)) {
    plot_height <- 6
    plot_width <- 8
    save_ggplot(p, output,
                width=plot_width,
                height=plot_height)
  }
  return(p)
}

data_dir <- argv$cov_dir

all_files <- list.files(data_dir)
filter_files <- all_files[grep("*.geneBodyCoverage.txt", all_files)]


split_str <- function(strings, Split) {
  for (i in 1:nchar(strings)) {
    if (substr(strings, i, i) == Split) {
      return(c(substr(strings, 1, i - 1), substr(strings, i + 1, nchar(strings))))
    }
  }
}


file_list <- list()
for (i in seq(length(filter_files))) {
  each_sample_df <- read.delim(paste(data_dir, filter_files[i], sep='/'))
  each_sample_df$reads_portion <- each_sample_df$count / sum(each_sample_df$count)
  each_sample_df$proportion <- prop.table(each_sample_df$count)
  sample_id <- split_str(filter_files[i], Split='.')[1]
  each_sample_df$sample <- sample_id
  each_sample_df[is.na(each_sample_df)] <- 0
  file_list[[i]] <- each_sample_df
  each_sample_out_name <- paste(sample_id, 'genebody_coverage', sep = '.')
  each_sample_out_path <- file.path(data_dir, each_sample_out_name)
  gene_cov_plot(each_sample_df, each_sample_out_path)
}

file_df <- ldply(file_list, data.frame)

samples <- unique(file_df$sample)
sample_number <- length(samples)
selected_num <- ifelse(sample_number < 9, sample_number, 9)
selected_df <- filter(file_df, sample %in% samples[1:selected_num])
plot_out <- file.path(data_dir, 'genebody_coverage.report')
selected_df$sample <- factor(selected_df$sample, levels = samples)
gene_cov_plot(selected_df, plot_out)
