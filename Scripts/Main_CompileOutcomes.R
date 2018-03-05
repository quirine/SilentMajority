library(mgcv)
library(rmutil)

# Abstract ----------------------------------------------------------------

# Our results indicate that 11 net infectiousness of people with asymptomatic infections is: 
test.1 <- c(AUC.asympto.1,AUC.asympto.2)
test.2 <- c(AUC.sympto.w.1,AUC.sympto.w.2)
median(test.1[order(test.1)]/test.2[order(test.2)])
HDIofMCMC(test.1[order(test.1)]/test.2[order(test.2)], 0.95)

# Due to their numerical prominence in the infectious reservoir, clinically inapparent infections in total could account for:
FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)]
FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI

# Of infections that ultimately result in symptoms, we estimate that 
XProps.brazil$Prop.sym.presym[length(foi.vec)]; XProps.brazil$Prop.sym.presym.low[length(foi.vec)]; XProps.brazil$Prop.sym.presym.up[length(foi.vec)]

# Only ... of DENV transmission is attributable to people with detected infections after they have developed symptoms. 
FOIProps.brazil$Prop.FOI.postsym.detected[length(foi.vec)]
FOIProps.brazil.CIs$Prop.FOI.postsym.detected.all.CI

# Author summary ----------------------------------------------------------

# At an individual level, we show that individuals with asymptomatic infections are capable of infecting ... as many mosquitoes as their symptomatic counterparts. 
median(test.1[order(test.1)]/test.2[order(test.2)])

# At a population level, we show that ..% of infections result from individuals who display no apparent symptoms at the time of transmission. 
(1-FOIProps.brazil$Prop.FOI.postsym[length(foi.vec)])


# Results -----------------------------------------------------------------

# The median net infectiousness to Ae. aegypti of A infections was lower than that of S infections, but of similar magnitude 
median(AUC.asympto.1[order(AUC.asympto.1)]/AUC.sympto.w.1[order(AUC.sympto.w.1)]); HDIofMCMC(AUC.asympto.1[order(AUC.asympto.1)]/AUC.sympto.w.1[order(AUC.sympto.w.1)])
median(AUC.asympto.2[order(AUC.asympto.2)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]); HDIofMCMC(AUC.asympto.2[order(AUC.asympto.2)]/AUC.sympto.w.2[order(AUC.sympto.w.2)])

# The median net infectiousness of 1° infections was greater than that of corresponding 2° infections 
median(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)])*100; HDIofMCMC(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)])*100
median(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]); HDIofMCMC(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)])

# Approximately one quarter of the net infectiousness of S infections occurred before symptom onset 
median(Prop.infectivity.prior.to.IIP.1); HDIofMCMC(Prop.infectivity.prior.to.IIP.1)
median(Prop.infectivity.prior.to.IIP.2); HDIofMCMC(Prop.infectivity.prior.to.IIP.2)
         
# we confirmed that A infections are more likely to be less infectious than S infections 
IntProbs[1,3]; IntProbs[2,4]
# and 2° infections are more likely to be less infectious than 1° infections 
IntProbs[2,1];IntProbs[4,3] 

# however, with both two-fold lesser or greater infectiousness compared to S infections appearing probable 
IntProbs.200[3,1]; IntProbs.200[4,2]
IntProbs.200[1,3]; IntProbs.200[2,4]

#1° A infections were not significantly less infectious than 2° S infections (Wilcoxon rank sum test, p = ...).
WSR.P$A.1vsS.2.p.value

# (A+IS):AS ratios between 1° and 2° infections do not differ significantly (likelihood ratio test, p = ...) and assessed these classes together:
#@@ pending: prim and sec
temp = summary(Post.exp.bb)
temp@coef[1,1]; temp@coef[1,1] - 1.96 * temp@coef[1,2]; temp@coef[1,1] + 1.96 * temp@coef[1,2]  

# Based on our metric of relative FoI, we estimated that ..% (CI: ..-..%) of human DENV infections are attributable to individuals 
#that do not present with apparent symptoms at the time when they are bitten by a susceptible mosquito 
(1-FOIProps.brazil$Prop.FOI.postsym[length(foi.vec)]); 
1 - FOIProps.brazil.CIs$Prop.FOI.postsym.CI

# We estimated that A and IS infections could be responsible for causing ..% (CI: ..-..%) of all human DENV infections
FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)]
FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI

# Of the remaining ..% (CI: ..-..%) of infections
1 - FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)]
1 - FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI

# ..% (CI: ..-..%) are attributable to bites by mosquitoes on people whose infection eventually becomes apparent, 
# and thus potentially detectable, after their onset of symptoms. 
XProps.brazil$Prop.sym.postsym[length(foi.vec)]; 
1 - XProps.brazil$Prop.sym.presym.low[length(foi.vec)]; 1 - XProps.brazil$Prop.sym.presym.up[length(foi.vec)]

# At a detection rate of 8% (24), ..% (CI: ..) of total infections were estimated to result from infected individuals after they are detected by surveillance systems
(1 - FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)]) * detection.rate
(1 - FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI) * detection.rate

#of total infections were estimated to result from infected individuals after they are detected by surveillance systems 
# (..% at detection rates of 5% and 15% (24), respectively).  
(1 - FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)]) * detection.rate.CI
(1 - FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI) * detection.rate.CI[1]
(1 - FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI) * detection.rate.CI[2]

# Sensitivity Analysis ----------------------------------------------------
# ). Under the assumption that the net infectiousness of post-secondary infections is equivalent to that of secondary infections, 
# we estimated the contribution of inapparent infections to be up to 
FOIProps.brazil.4$Prop.FOI.inap.asym[length(foi.vec)] / FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)] - 1
FOIProps.brazil.4.CIs$Prop.FOI.inap.asym.CI / FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI  - 1
# higher than if post-secondary infections had not been assumed to contribute to transmissiozn (Fig S1).

# Under the assumption that IS infections are more similar in their infectiousness to A than to AS infections, the estimated contribution 
# of inapparent infections was reduced from ... 
FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)]
# to ..., reflecting a lower bound on this assumption.
X.temp <- Disease.prop.calculator(infec.hist.brazil.small,c(inapp.ratio.rev[1],1-inapp.ratio.rev[1],0), c(inapp.ratio.rev[2],1-inapp.ratio.rev[2],0), c(inapp.ratio.rev[3],1-inapp.ratio.rev[3],0)) 
temp <- Sympto.FOI.sampler(prev = infec.hist.brazil$prev.vec2,X = X.temp[1,], 2, SA.foi.vec[1], Ca.all, Cm.all, Cs.all, cores = 3)
1 - median(temp)
1 - round( HDIofMCMC(temp) , 2)

# The impact of accounting for the differential viral trajectories of severe dengue cases (18) was minor due to their small numerical prominence 
# (28), but their inclusion did increase the contribution of post-symptomatic DAS infections from 
FOIProps.brazil$Prop.FOI.postsym.detected[length(foi.vec)]
FOIProps.brazil.CIs$Prop.FOI.postsym.detected.all.CI

FOIProps.brazil.DHF$Prop.FOI.postsym.detected[length(foi.vec)]
FOIProps.brazil.DHF.CIs$Prop.FOI.postsym.detected.all.CI

# Although variability in As:IS:AS ratios likely does drive some variability in the contribution of As+IS infections, 
# our analysis suggests that their contribution is unlikely to be less than ..% of their numerical prominence (Fig 5)  
1 - max(df$Res.inap.as)

# Methods -----------------------------------------------------------------
# As:(IS+AS) ratio CI:
0.092 + qnorm(1-0.05/2)*sqrt(0.092*(1-0.092)/141)
0.092 - qnorm(1-0.05/2)*sqrt(0.092*(1-0.092)/141)






# Archive -----------------------------------------------------------------



paste(c('the proportion of infections attributable to inap and asym at FOI= 0.1 is',FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)] ))
paste(c('with 95%, median, and 5% bounds of', FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI))

paste(c(FOIProps.brazil$Prop.FOI.asym[length(foi.vec)], '% of which by asymptomatic infections'))
paste(c('with 95%, median, and 5% bounds of', FOIProps.brazil.CIs$Prop.FOI.asym.CI))

paste(c('About ',XProps.brazil$Prop.sym.presym[length(foi.vec)],'% of the remaining infections are estimated to result from individuals before they present with clinical symptoms.'   ) )
paste(c('with 95%, median, and 5% bounds of', XProps.brazil$Prop.sym.presym.low[length(foi.vec)], XProps.brazil$Prop.sym.presym.up[length(foi.vec)]) )

paste(c('This equates to ',(1-FOIProps.brazil$Prop.FOI.postsym[length(foi.vec)]), '% of infections annually being attributable to infected individuals that have not (yet) presented with apparent symptoms. '))
paste(c('with 95%, median, and 5% bounds of', (1-FOIProps.brazil.CIs$Prop.FOI.postsym.CI) ))


# Main text
# viremia adjustment factor
quantile(Ratios.asym, c(0.25,0.5,0.75))

# disease proportions
paste(c('We find that Asym:Sym ratios vary substantially across regions (SFig1a), 
    and that secondary infections are more likely to result in symptomatic infection 
     (',asym.ratio[2]$pred ,' % in primary vs.', asym.ratio[1]$pred,'% in secondary).' ))

paste( c('with CI bounds of',asym.ratio[2]$ci.lb, 'and' , asym.ratio[2]$ci.ub ))
paste( c('and',asym.ratio[1]$ci.lb, 'and' , asym.ratio[1]$ci.ub ))

paste( c('The symptomatic infections can further be divided into inapparent and 
   apparent using the Inap:Ap-ratio (', inapp.ratio$pred,'(',inapp.ratio$ci.lb,'-',inapp.ratio$ci.ub,') % vs.',
         1-inapp.ratio$pred,'(',1-inapp.ratio$ci.lb,'-',1-inapp.ratio$ci.ub,')%,pre-exposure history insignificant
         (p=',res.exp$pval[2] ,'), SFig1b)')) 

paste( c('We estimate that the proportion of infections to be asymptomatic varies between', Xa.norm[length(foi.vec)],'-',Xa.norm[1],'%, depending on the prevailing FOI (Figure 1).      '))

# viremia and infectivity
# paste( c('the overall infectivity of individuals with primary DENV is twice as high (',mean(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), HPDinterval(mcmc(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), 0.95) ),
#          '%) compared to secondary infections for symptomatic infections, and three times as high for asymptomatic infections (',
# mean(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)]), HPDinterval(mcmc(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)]), 0.95),'%')

paste( c('On average asymptomatic individuals are estimated to infect mosquitoes at a lower probability (',median(AUC.asympto.1[order(AUC.asympto.1)]/AUC.sympto.w.1[order(AUC.sympto.w.1)]), quantile(AUC.asympto.1[order(AUC.asympto.1)]/AUC.sympto.w.1[order(AUC.sympto.w.1)], c(0.25,0.75)),'
         % (CI) for primary and ',median(AUC.asympto.2[order(AUC.asympto.2)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), quantile(AUC.asympto.2[order(AUC.asympto.2)]/AUC.sympto.w.2[order(AUC.sympto.w.2)], c(0.25,0.75)),'% (CI) secondary infections) than mild and severe symptomatic infections with the same pre-exposure history (Figure 2). '))

paste(c('Among A infections, the net infectiousness of primary infections is double that of secondary infections (',median(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)]), quantile(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)], c(0.25,0.75)) ,')'))

paste(c('The net infectiousness of primary S infections is triple that of secondary infections ',median(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), quantile(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)], c(0.25,0.75)) ,')'))

paste( c('The largest proportion of the overall infectivity of symptomatic infections (',median(Prop.infectivity.prior.to.IIP.1),quantile(Prop.infectivity.prior.to.IIP.1,c(0.25,0.75))
         ,'% (CI) for primary and', median(Prop.infectivity.prior.to.IIP.2),quantile(Prop.infectivity.prior.to.IIP.2,c(0.25,0.75)),
         '% (CI) for secondary) occurs before symptom onset.'))

d <- hist((AUC.asympto.1/ AUC.sympto.w.1[order(AUC.sympto.w.1)]), breaks = 100, plot = FALSE) 
d$density <- d$density*diff(d$breaks)[1]
paste(c('of primary A infections, there is a probability of', sum(d$density[1:5]),'% that they contribute less than 50% relative to I-prim'))
paste(c('and', sum(d$density[-seq(0,20)]),'% that they are more than twice as infectious'))
d <- hist((AUC.asympto.2/ AUC.sympto.w.2[order(AUC.sympto.w.2)]), breaks = 100, plot = FALSE) 
d$density <- d$density*diff(d$breaks)[1]
paste(c('of secondary A infections, there is a probability of', sum(d$density[1]),'% that they contribute less than 50% relative to I-prim'))
paste(c('and', sum(d$density[-seq(0,4)]),'% that they are more than twice as infectious'))

paste( c(' Primary A infections may well be more infectious than secondary S infection:', A.1vsS.2[c(1,3)]))

# FOI contributions
paste(c('This equates to ',(1-FOIProps.brazil$Prop.FOI.postsym[length(foi.vec)]), '% of infections annually being attributable to infected individuals that have not (yet) presented with apparent symptoms. '))
paste(c('with 95%, median, and 5% bounds of', (1-FOIProps.brazil.CIs$Prop.FOI.postsym.CI) ) )

paste( c('we find asymptomatic infections to be responsible for a median of', FOIProps.brazil$Prop.FOI.asym[length(foi.vec)],
         '% (',FOIProps.brazil.CIs$lower.1, FOIProps.brazil.CIs$upper.1,') of human DENV infections and inapparent symptomatic infections', FOIProps.brazil$Prop.FOI.inap[length(foi.vec)],
         '% (',FOIProps.brazil.CIs$lower.2, FOIProps.brazil.CIs$upper.2,'(Figure 3)'))
paste( c('These estimates are found to be relatively robust across transmission settings with their means ranging from ', FOIProps.brazil$Prop.FOI.asym[1], FOIProps.brazil$Prop.FOI.asym[length(foi.vec)] ))

paste( c('we find inapparent and asymptomatic infections to be responsible for on average', FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)],
         '% (',FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI,') (Figure 3)'))
paste( c('These estimates are found to be relatively robust across transmission settings with their means ranging from ', FOIProps.brazil$Prop.FOI.inap.asym[1], FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)] ))


paste( c('Of the remaining',1-FOIProps.brazil$Prop.FOI.asym[length(foi.vec)]-FOIProps.brazil$Prop.FOI.inap[length(foi.vec)],'(',1 - FOIProps.brazil.CIs$Prop.FOI.inap.asym.CI,
         ')% of the infections, ',XProps.brazil$Prop.sym.postsym[length(foi.vec)],'(',1-XProps.brazil$Prop.sym.presym[length(foi.vec)],1-XProps.brazil$Prop.sym.presym.low[length(foi.vec)], 1-XProps.brazil$Prop.sym.presym.up[length(foi.vec)] ,
         ')% of the transmission events are likely to occur after symptoms have presented, 
         equating to just ',FOIProps.brazil$Prop.FOI.postsym[length(foi.vec)],'(',FOIProps.brazil.CIs$lower.5, FOIProps.brazil.CIs$upper.5,') % of total infections.'))
paste( c('We find this proportion to be even higher in settings with a population fully na??ve to DENV  (',FOIProps.brazil.emerging$Prop.FOI.asym[length(foi.vec)]
         ,'(',FOIProps.brazil.emerging.CIs$lower.1, FOIProps.brazil.emerging.CIs$upper.1,')%) '))
paste( c('We find this proportion to be even higher in settings with a population fully na??ve to DENV  (',FOIProps.brazil.emerging$Prop.FOI.inap.asym[length(foi.vec)]
         ,'(',FOIProps.brazil.emerging.CIs$lower.3, FOIProps.brazil.emerging.CIs$upper.3,')%) '))

# detected

paste('Only ... of total infections results from detected individuals after they present at a health clinic')
FOIProps.brazil$Prop.FOI.postsym.detected[length(foi.vec)]
FOIProps.brazil.CIs$Prop.FOI.postsym.detected.all.CI
# + CIs from Stanaway
FOIProps.brazil$Prop.FOI.postsym.detected.low[length(foi.vec)]
FOIProps.brazil.CIs$Prop.FOI.postsym.detected.low.all.CI

FOIProps.brazil$Prop.FOI.postsym.detected.high[length(foi.vec)]
FOIProps.brazil.CIs$Prop.FOI.postsym.detected.high.all.CI

# sensitivity analysis
idx <- which.min(abs(Prop.FOI.All.brazil.inap[length(foi.vec),] - median(Prop.FOI.All.brazil.inap[length(foi.vec),])))
paste( c('assuming inapparents to be equally infectious as asymptomatics, lowers there overall contribution to',
         1-Prop.FOI.All.brazil.inap[length(foi.vec),idx], '(', quantile(1-Prop.FOI.All.brazil.inap[length(foi.vec),],c(0.25,0.75)),')' ))

paste( c('We demonstrate that the contribution of asymptomatic infections is up to ',FOIProps.brazil.4$Prop.FOI.inap.asym[length(foi.vec)] - FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)],
         '(',quantile(FOIProps.brazil.4$Prop.FOI.inap.asym - FOIProps.brazil$Prop.FOI.inap.asym, c(0.25,0.75)),')% higher when accounting for these infections'))


# archive - Cred Intervals ------------------------------------------------
# 
# # Abstract
# 
# paste(c('the proportion of infections attributable to inap and asym at FOI= 0.1 is',FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)] ))
# paste(c('with 95%, median, and 5% bounds of', FOIProps.brazil.CIs$lower.3, FOIProps.brazil.CIs$upper.3))
# 
# paste(c(FOIProps.brazil$Prop.FOI.asym[length(foi.vec)], '% of which by asymptomatic infections'))
# paste(c('with 95%, median, and 5% bounds of', FOIProps.brazil.CIs$lower.1, FOIProps.brazil.CIs$upper.1))
# 
# paste(c('About ',XProps.brazil$Prop.sym.presym[length(foi.vec)],'% of the remaining infections are estimated to result from individuals before they present with clinical symptoms.'   ) )
# paste(c('with 95%, median, and 5% bounds of', XProps.brazil$Prop.sym.presym.low[length(foi.vec)], XProps.brazil$Prop.sym.presym.up[length(foi.vec)]) )
# 
# paste(c('This equates to ',(1-FOIProps.brazil$Prop.FOI.postsym[length(foi.vec)]), '% of infections annually being attributable to infected individuals that have not (yet) presented with apparent symptoms. '))
# paste(c('with 95%, median, and 5% bounds of', (1-FOIProps.brazil.CIs$lower.5), (1-FOIProps.brazil.CIs$upper.5) ) )
# 
# # Main text
# # disease proportions
# paste(c('We find that Asym:Sym ratios vary substantially across regions (SFig1a), 
#         and that secondary infections are more likely to result in symptomatic infection 
#         (',asym.ratio[2]$pred ,' % in primary vs.', asym.ratio[1]$pred,'% in secondary).' ))
# 
# paste( c('with CI bounds of',asym.ratio[2]$ci.lb, 'and' , asym.ratio[2]$ci.ub ))
# paste( c('and',asym.ratio[1]$ci.lb, 'and' , asym.ratio[1]$ci.ub ))
# 
# paste( c('The symptomatic infections can further be divided into inapparent and 
#          apparent using the Inap:Ap-ratio (', inapp.ratio$pred,'(',inapp.ratio$ci.lb,'-',inapp.ratio$ci.ub,') % vs.',
#          1-inapp.ratio$pred,'(',1-inapp.ratio$ci.lb,'-',1-inapp.ratio$ci.ub,')%,pre-exposure history insignificant
#          (p=',res.exp$pval[2] ,'), SFig1b)')) 
# 
# paste( c('We estimate that the proportion of infections to be asymptomatic varies between', Xa.norm[length(foi.vec)],'-',Xa.norm[1],'%, depending on the prevailing FOI (Figure 1).      '))
# 
# # viremia and infectivity
# # paste( c('the overall infectivity of individuals with primary DENV is twice as high (',mean(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), HPDinterval(mcmc(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), 0.95) ),
# #          '%) compared to secondary infections for symptomatic infections, and three times as high for asymptomatic infections (',
# # mean(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)]), HPDinterval(mcmc(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)]), 0.95),'%')
# 
# paste( c('On average asymptomatic individuals are estimated to infect mosquitoes at a lower probability (',mean(AUC.asympto.1[order(AUC.asympto.1)]/AUC.sympto.w.1[order(AUC.sympto.w.1)]), HPDinterval(mcmc(AUC.asympto.1[order(AUC.asympto.1)]/AUC.sympto.w.1[order(AUC.sympto.w.1)]), 0.95),'
#          % (CI) for primary and ',mean(AUC.asympto.2[order(AUC.asympto.2)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), HPDinterval(mcmc(AUC.asympto.2[order(AUC.asympto.2)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), 0.95),'% (CI) secondary infections) than mild and severe symptomatic infections with the same pre-exposure history (Figure 2). '))
# 
# paste(c('Among A infections, the net infectiousness of primary infections is double that of secondary infections (',mean(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)]), HPDinterval(mcmc(AUC.asympto.1[order(AUC.asympto.1)]/AUC.asympto.2[order(AUC.asympto.2)]), 0.95) ,')'))
# 
# paste(c('The net infectiousness of primary S infections is triple that of secondary infections ',mean(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), HPDinterval(mcmc(AUC.sympto.w.1[order(AUC.sympto.w.1)]/AUC.sympto.w.2[order(AUC.sympto.w.2)]), 0.95) ,')'))
# 
# paste( c('The largest proportion of the overall infectivity of symptomatic infections (',mean(Prop.infectivity.prior.to.IIP.1),HPDinterval(mcmc(Prop.infectivity.prior.to.IIP.1),0.95)
#          ,'% (CI) for primary and', mean(Prop.infectivity.prior.to.IIP.2),HPDinterval(mcmc(Prop.infectivity.prior.to.IIP.2),0.95),
#          '% (CI) for secondary) occurs before symptom onset.'))
# 
# paste( c(' Primary A infections may well be more infectious than secondary S infection:', A.1vsS.2[c(1,3)]))
# 
# # FOI contributions
# paste(c('This equates to ',(1-FOIProps.brazil$Prop.FOI.postsym[length(foi.vec)]), '% of infections annually being attributable to infected individuals that have not (yet) presented with apparent symptoms. '))
# paste(c('with 95%, median, and 5% bounds of', (1-FOIProps.brazil.CIs$lower.5), (1-FOIProps.brazil.CIs$upper.5) ) )
# 
# paste( c('we find asymptomatic infections to be responsible for on average', FOIProps.brazil$Prop.FOI.asym[length(foi.vec)],
#          '% (',FOIProps.brazil.CIs$lower.1, FOIProps.brazil.CIs$upper.1,') of human DENV infections and inapparent symptomatic infections', FOIProps.brazil$Prop.FOI.inap[length(foi.vec)],
#          '% (',FOIProps.brazil.CIs$lower.2, FOIProps.brazil.CIs$upper.2,'(Figure 3)'))
# paste( c('These estimates are found to be relatively robust across transmission settings with their means ranging from ', FOIProps.brazil$Prop.FOI.asym[1], FOIProps.brazil$Prop.FOI.asym[length(foi.vec)] ))
# 
# paste( c('we find inapparent and asymptomatic infections to be responsible for on average', FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)],
#          '% (',FOIProps.brazil.CIs$lower.3, FOIProps.brazil.CIs$upper.3,') (Figure 3)'))
# paste( c('These estimates are found to be relatively robust across transmission settings with their means ranging from ', FOIProps.brazil$Prop.FOI.inap.asym[1], FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)] ))
# 
# 
# paste( c('Of the remaining',1-FOIProps.brazil$Prop.FOI.asym[length(foi.vec)]-FOIProps.brazil$Prop.FOI.inap[length(foi.vec)],'(',1 - FOIProps.brazil.CIs$lower.3,1- FOIProps.brazil.CIs$upper.3 ,
#          ')% of the infections, ',XProps.brazil$Prop.sym.postsym[length(foi.vec)],'(',1-XProps.brazil$Prop.sym.presym[length(foi.vec)],1-XProps.brazil$Prop.sym.presym.low[length(foi.vec)], 1-XProps.brazil$Prop.sym.presym.up[length(foi.vec)] ,
#          ')% of the transmission events are likely to occur after symptoms have presented, 
#          equating to just ',FOIProps.brazil$Prop.FOI.postsym[length(foi.vec)],'(',FOIProps.brazil.CIs$lower.5, FOIProps.brazil.CIs$upper.5,') % of total infections.'))
# paste( c('We find this proportion to be even higher in settings with a population fully na??ve to DENV  (',FOIProps.brazil.emerging$Prop.FOI.asym[length(foi.vec)]
#          ,'(',FOIProps.brazil.emerging.CIs$lower.1, FOIProps.brazil.emerging.CIs$upper.1,')%) '))
# paste( c('We find this proportion to be even higher in settings with a population fully na??ve to DENV  (',FOIProps.brazil.emerging$Prop.FOI.inap.asym[length(foi.vec)]
#          ,'(',FOIProps.brazil.emerging.CIs$lower.3, FOIProps.brazil.emerging.CIs$upper.3,')%) '))
# 
# 
# # sensitivity analysis
# paste( c('assuming inapparents to be equally infectious as asymptomatics, lowers there overall contribution to',
#          1-mean(Prop.FOI.All.brazil.inap[length(foi.vec),]), '(', HPDinterval(mcmc(1-Prop.FOI.All.brazil.inap[length(foi.vec),]),0.95),')' ,)
#        
#        paste( c('We demonstrate that the contribution of asymptomatic infections is up to ',FOIProps.brazil.4$Prop.FOI.inap.asym[length(foi.vec)] - FOIProps.brazil$Prop.FOI.inap.asym[length(foi.vec)],
#                 '(',HPDinterval(mcmc(FOIProps.brazil.4$Prop.FOI.inap.asym - FOIProps.brazil$Prop.FOI.inap.asym), 0.95),')% higher when accounting for these infections'))
#        
