suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(omplotr))
suppressMessages(library(argparser))


p <- arg_parser('rmats analysis summary')
p <- add_argument(p,'--as_dir',help = 'rmats result directory.')
argv <- parse_args(parser = p)

as_types <- c('SE', 'A5SS', 'A3SS', 'MXE', 'RI')
as_dir <- argv$as_dir
compare_name <- basename(as_dir)

filter_data <- function(as_type, as_dir=NULL, diff_cut=0.05, reads_type='JC') {
  as_file_name <- paste(as_type, 'MATS', reads_type,'txt', sep='.')
  as_file <- file.path(as_dir, as_file_name)
  as_df = fread(as_file)
  as_df$as_type <- as_type
  as_df$reads_type <- reads_type
  f_as_df <- as_df[as_df$FDR <= diff_cut]
  f_as_df$direction <- ifelse(f_as_df$IncLevelDifference > 0, 'up', 'down')
  out_df <- f_as_df[, c('direction', 'as_type', 'GeneID', 'chr', 'reads_type')]
  return(out_df)
}

as_df_list <- lapply(as_types, filter_data, as_dir=as_dir)
merged_as_df <- ldply(as_df_list, data.frame)
jcec_as_df_list <- lapply(as_types, filter_data, as_dir=as_dir, reads_type='JCEC')
jcec_merged_as_df <- ldply(jcec_as_df_list, data.frame)
all_merged_as_df <- rbind(merged_as_df, jcec_merged_as_df)
as_count <- data.frame(table(all_merged_as_df$as_type, 
                             all_merged_as_df$direction,
                             all_merged_as_df$reads_type))

p <- ggplot(as_count, aes(Var1, Freq, fill = Var2)) + 
  geom_bar(position = 'dodge', color = 'black', stat = 'identity') +
  theme_onmath() +
  scale_fill_brewer(palette = 'Set1') +
  xlab('') + ggtitle(compare_name) +
  facet_wrap(~Var3) +
  guides(fill=guide_legend(title = '')) +
  ylab('AS Event Number') +
  xlab('')

out_prefix <- file.path(as_dir, paste(compare_name, 'as_summary', sep = '.'))
save_ggplot(ggplot_out = p,output = out_prefix,
            width = 8, height = 6)
