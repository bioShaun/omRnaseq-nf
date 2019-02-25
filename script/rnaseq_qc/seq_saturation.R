options(warn = -1)
options(stringsAsFactors = F)
suppressMessages(library('reshape2'))
suppressMessages(library('argparser'))
suppressMessages(library('plyr'))
suppressMessages(library('dplyr'))
suppressMessages(library('omplotr'))

p <- arg_parser('seq saturation')
p <- add_argument(p,'--saturation_count',help = 'saturation count directory.')
argv <- parse_args(parser = p)


seq_saturation <- function(filter_df, output) {
  sample_number <- length(unique(filter_df$sample))
  
  p <- ggplot(filter_df, aes(reads_portion, trancripts, color=sample, group=sample)) +
    geom_point() +
    geom_line() +
    theme_onmath() +
    scale_color_brewer(palette = 'Set1') +
    xlab('Proportion of Reads') + ylab('Transcripts at 10x')
  p
  
  if (sample_number > 1) {
    facet_wrap_ncol = round(sqrt(sample_number))
    p <- p + guides(color=guide_legend(title = ''))
  } else {
    p <- p + guides(color=F)
  }
  
  if (! is.null(output)) {
    plot_height <- 6
    plot_width <- 12
    save_ggplot(p, output,
                width=plot_width,
                height=plot_height)
  }
  return(p)
}

data_dir <- argv$saturation_count

all_files <- list.files(data_dir)
filter_files <- all_files[grep("*.rawCount.xls", all_files)]


split_str <- function(strings, Split) {
  for (i in 1:nchar(strings)) {
    if (substr(strings, i, i) == Split) {
      return(c(substr(strings, 1, i - 1), substr(strings, i + 1, nchar(strings))))
    }
  }
}


file_list <- list()
for (i in seq(length(filter_files))) {
  each_sample_df <- read.delim(paste(data_dir, filter_files[i], sep='/'), check.names = F)
  count_filter_df <- t(each_sample_df[, -seq(1:6)] > 10)
  passed_cut_num <- as.data.frame(rowSums(count_filter_df))
  colnames(passed_cut_num) <- c('trancripts')
  sample_id <- split_str(filter_files[i], Split='.')[1]
  passed_cut_num$sample <- sample_id
  passed_cut_num$reads_portion <- rownames(passed_cut_num)
  
  passed_cut_num$reads_portion <- factor(passed_cut_num$reads_portion,
                                         levels = passed_cut_num$reads_portion)
  file_list[[i]] <- passed_cut_num
  each_sample_out_name <- paste(sample_id, 'sequence_saturation', sep = '.')
  each_sample_out_path <- file.path(data_dir, each_sample_out_name)
  seq_saturation(passed_cut_num, each_sample_out_path)
}

file_df <- ldply(file_list, data.frame)

samples <- unique(file_df$sample)
sample_number <- length(samples)
selected_num <- ifelse(sample_number < 9, sample_number, 9)
selected_df <- filter(file_df, sample %in% samples[1:selected_num])
plot_out <- file.path(data_dir, 'sequence_saturation.report')
selected_df$sample <- factor(selected_df$sample, levels = samples)
seq_saturation(selected_df, plot_out)
