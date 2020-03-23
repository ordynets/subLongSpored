# Load required packages
library(here)
library(dplyr)

# Set the current directory as working directory
set_here()

# Load spore table and adjust it
sub <- read.csv ("SUB_spore_measurements_20191011_v1.csv",  sep=",")

# Establish lower and upper quantiles
prob=c(0.05, 0.95)

spores.specim.summ <- sub %>% 
  group_by(Specimen_ID) %>% 
  summarise(Lmin=min(L_sp),
            L05=quantile(L_sp, prob=0.05),
            Lmean=mean(L_sp),
            L95=quantile(L_sp, prob=0.95),
            Lmax=max(L_sp),
            
            Wmin=min(W_sp),
            W05=quantile(W_sp, prob=0.05),
            Wmean=mean(W_sp),
            W95=quantile(W_sp, prob=0.95),
            Wmax=max(W_sp),
            
            Qmin=min(Q_sp),
            Q05=quantile(Q_sp, prob=0.05),
            Qmean=mean(Q_sp),
            Q95=quantile(Q_sp, prob=0.95),
            Qmax=max(Q_sp),
            
            Spore_count=n())


# Write result in a excel file
write.table(spores.specim.summ, file="Spore_ranges_specimens_20191011_v0.csv",sep=",",quote=F)
