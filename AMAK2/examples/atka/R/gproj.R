getwd()
bdf <- read.table("../proj/BSAI_atka_out/bigfile.out",header=T)
bdft <- data.frame(bdf) %>% filter(Alternative==1) %>% slice(1:420) %>% transmute(Yr,SSB,sim=sort(rep(1:30,14)))
df <- read.table("../proj/BSAI_atka_out/percentdb.out",header=F,as.is=T)
names(df) <- c("mod","alt","yr","idx","val")
tbl_df(bdft)
tbl_df(tdf)
tail(bdft)
## SSB plot...
filter(df,grepl('SSB', idx),alt==1) %>% ggplot(aes(x=yr,y=val,colour=as.factor(idx))) + 
     geom_line(size=2) + mytheme + ylab("Spawning biomass") + xlab("Year") + ylim(c(0,350))
tbl_df(tdf)
tdf <- filter(df,grepl('SSBFabc|SSBFofl|SSBLCI|SSBUCI|SSBMedi', idx),alt==1) %>% mutate(C=as.factor(idx),C=factor(C,levels=c("SSBFabc","SSBFofl","SSBMedian","SSBLCI","SSBUCI")))
ggplot(tdf,aes(x=yr,y=val,col=C, linetype=C)) + 
    geom_line(size=1.4) + mytheme + ylab("Spawning biomass (kt)") + xlab("Year") + ylim(c(0,200)) +
		scale_size_manual(values = c(rep(1.3,3),rep(0.8,2))) +
		scale_linetype_manual(values = c(rep("solid",3),rep("dotted",2),rep("solid",40))) +
		scale_color_manual(   values = c("green","red","blue",rep("grey",42)))  + theme(legend.position="none")  +
		geom_line(data=bdft,aes(x=Yr,y=SSB,group=sim))
		,colour="grey")
		bdft$SSB
		bdft$sim

## Catch plot...
tdf <- filter(df,grepl('Cabc|CLC|Cofl|CUCI|CMedi', idx),alt==1) %>% mutate(C=as.factor(idx),C=factor(C,levels=c("Cabc","Cofl","CMedian","CLCI","CUCI")))
ggplot(tdf,aes(x=yr,y=val,col=C, linetype=C)) + 
    geom_line(size=1.4) + mytheme + ylab("Catch (kt)") + xlab("Year") + ylim(c(0,200)) +
		scale_size_manual(values = c(rep(1.3,3),rep(0.8,2))) +
		scale_linetype_manual(values = c(rep("solid",3),rep("dotted",2))) +
		scale_color_manual(   values = c("green","red","blue","grey","grey"))  +theme(legend.position="none")