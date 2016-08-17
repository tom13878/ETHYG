# OBSOLETE CODE

# Note that MPP cannot be calculated for plots with N=0 and are therefore set to 0.
db2_olsCD1 <- db1 %>% mutate(elastfert = ifelse(AEZ2=="Tropic-warm", coef(model)["AEZ2Tropic-warm:logN"], 
                                                ifelse(AEZ2=="Tropic-cool / semi-arid", coef(model)["AEZ2Tropic-cool / semi-arid:logN"],
                                                       ifelse(AEZ2=="Tropic-cool / sub-humid", coef(model)["AEZ2Tropic-cool / sub-humid:logN"], 
                                                              coef(model)["AEZ2Tropic-cool:logN"]))),
                             phconstant = ifelse(phdum==2, coef(model)["phdum22"], 
                                                 ifelse(phdum==3, coef(model)["phdum23"],0)),
                             lnA = coef(model)["(Intercept)"] +
                               (coef(model)["noN"]*noN) +
                               (coef(model)["logarea"]*logarea) +
                               (coef(model)["hybrd"]*hybrd) +
                               (coef(model)["manure"]*manure) +
                               (coef(model)["herb"]*herb) +
                               (coef(model)["fung"]*fung) +
                               (coef(model)["legume"]*legume) +
                               (coef(model)["irrig"]*irrig) +
                               (coef(model)["slope"]*slope) +
                               (coef(model)["elevation"]*elevation) +
                               (coef(model)["SOC2"]*SOC2) +
                               (phconstant) +
                               (coef(model)["rain_wq"]* rain_wq) +
                               (coef(model)["rain_wq_2"]*rain_wq_2) +
                               (coef(model)["crop_count2"]*crop_count2),
                             lnA2 = lnA,
                             constantfactor = exp(lnA2)*elastfert*(lab^coef(model)["loglab"]),
                             MPP= ifelse(N==0,NA,exp(lnA2)*elastfert*(lab^coef(model)["loglab"])*(N^(elastfert-1))),
                             Npm = (Pn/(constantfactor*Pm))^(1/(elastfert-1)),
                             Ndif = N-Npm)                       


sumzone_olsCD1<- db2_olsCD1 %>% group_by(region_lsms) %>%
  summarize(
    Ncon=mean(ifelse(N>0, N, NA), na.rm = T),
    N=mean(N, na.rm=T),
    Npm=mean(Npm, na.rm=T),
    MPPmean=mean(MPP[!is.infinite(MPP)], na.rm=T),
    MVC=mean((Pm[!is.infinite(MPP)]*MPP[!is.infinite(MPP)])/Pn[!is.infinite(MPP)], na.rm=T),
    Ndif=mean(Ndif, na.rm=T),
    Number=n())
