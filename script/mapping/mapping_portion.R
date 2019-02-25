suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(omplotr))
suppressMessages(library('argparser'))
options(stringsAsFactors = F)


p <- arg_parser('Mapping rate plot')
p <- add_argument(p,'--mapping_file',help = 'Mapping overall stats file.')
p <- add_argument(p,'--out_prefix',help = 'Plot output prefix.')
argv <- parse_args(parser = p)

mapping_file <- argv$mapping_file
out_prefix <- argv$out_prefix

mapping_df <- read.delim(mapping_file)

geome_portion_df <- mapping_df[, seq(3,5)] * mapping_df$Mapping_rate
geome_portion_df$sample_id <- mapping_df$Sample_id
m_geome_portion_df <- melt(geome_portion_df, id.vars = 'sample_id')
m_geome_portion_df$variable <- factor(m_geome_portion_df$variable,
                                      levels = c('Intergenic',
                                                 'Intron',
                                                 'Exon'))

sample_num <- dim(mapping_df)[1]

plot_width <- 6 + sample_num * 0.2
plot_height <- 6 + sample_num * 0.1

p <- ggplot(m_geome_portion_df, aes(sample_id, value, fill=variable)) +
  geom_bar(stat = 'identity', color='white') +
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0,1,0.1)) +
  scale_fill_brewer(palette = 'Set1') +
  theme_onmath() + theme(axis.text.x = element_text(angle = 45,
                                                    hjust = 1)) +
  guides(fill = guide_legend(title = '')) +
  xlab('') + ylab('Mapping rate')

save_ggplot(p, out_prefix,
            width = plot_width,
            height = plot_height)
