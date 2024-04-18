
#---Multiple SS models------------------
ssmods <- SSgetoutput(dirvec = c("ss/base","ss/mod","ss/mix","ss/high"))

r1 <- SS_obj(SS_output(dir = here("ss","base")),src="base") 
r2 <- SS_obj(SS_output(dir = here("ss","mod")),src="mod")  
r3 <- SS_obj(SS_output(dir = here("ss","mix")),src="mix")  
r4 <- SS_obj(SS_output(dir = here("ss","high")),src="high")  
#r5 <- SS_obj(SS_output(dir = here("ss","autocor")),src="autocor")  
ss_sel <- rbind(r1$sel,r2$sel,r3$sel,r4$sel)#,r5$sel)
(r1$sel)
myl<-list(c(r1,r2,r3,r4))
#Q: use lapply to apply function to all in my list of models 
fun <- function(x) {
  sel <- x$sel
  compute_matrix_summary(sel[,2:16])
}
fun(r1)

#

lapply(myl,fun)
compute_matrix_summary(r1$sel[,2:16])
compute_matrix_summary(r2$sel[,2:16])
Plot_Sel(sel=ss_sel,lage=15)
Plot_Sel_age(sel=ss_sel,lage=15)
all_ts <- rbind( sam_run$ts,pm_run$ts,am_run$ts,ss_run$ts,gp_run$ts)
ss_ts  <- rbind(r1$ts,r2$ts,r3$ts,r4$ts)#,r5$sel)
Plot_SSB(ss_ts)
?SSgetoutput
summary1 <- SSsummarize(ssmods)
names(summary1)
summary1$modelnames<- c("base","mod","mix","high")
SSplotComparisons(summary1)
print(summary1)

gt::gt(summary1$likelihoods)

