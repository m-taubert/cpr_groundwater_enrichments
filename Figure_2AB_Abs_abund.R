# Figure 2A and 2B

library(ggplot2)
# Figure 2A
# required data
F2AB.raw <- meta[which(meta$filter.fraction==0.2),c("Csource","Period","qPCR","CPR.rel.abund")] # CPR.rel.abund is added in script for Fig 1A, run it before this script

# add in situ data
F2AB <- rbind(F2AB.raw, in.situ[,c("Csource","Period","qPCR","CPR.rel.abund")])

# calculate absolute abundances for Fig 2B
F2AB$CPR.abs.abund <- F2AB$qPCR*F2AB$CPR.rel.abund

# plot Figure 2A
ggplot(F2AB, aes(x=Period, y=log10(qPCR), fill=Csource)) + 
 geom_boxplot(position=position_dodge(width=0.8), outlier.colour="black", outlier.size=3) +
 labs(x = "Period", y = "Abs. abundance / gene copies L-1", fill = "Treatment"
 ) +
 scale_fill_manual(values = c("#B69E8B", "#808080", "#3794bf", "#df8640"), 
                   name = "",
                   breaks = c("in situ", "none", "autotrophy", "methylotrophy"),
                   labels = c("in situ", "none", "autotrophy", "methylotrophy")) +
 scale_x_discrete(limits=c("in situ", "day0","1week","2weeks","3weeks",">3weeks")) +
 #coord_cartesian(ylim=c(4,12)) +
 scale_y_continuous(limits = c(4, 11), breaks = seq(4, 10, by = 2)) +
 theme(legend.position="bottom")

# plot Figure 2B

ggplot(F2AB, aes(x = Period, y = log10(CPR.abs.abund), fill = Csource)) + 
 geom_boxplot() +
 scale_fill_manual(values = c("#B69E8B","#808080", "#3794bf", "#df8640"), 
                   breaks = c("in situ", "none", "autotrophy", "methylotrophy"),
                   labels = c("in situ", "None", "Autotrophy", "Methylotrophy")) +
 scale_x_discrete(limits=c("in situ", "day0","1week","2weeks","3weeks",">3weeks")) +
 scale_y_continuous(limits = c(4, 11), breaks = seq(4, 10, by = 2)) +
 labs(y = "Absolute Abundance of CPR \n gene copies / 1 liter", x = "Incubation Period", fill = "Treatment")+
  theme(legend.position = "bottom")

