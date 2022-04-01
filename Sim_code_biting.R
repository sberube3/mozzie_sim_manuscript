library(dplyr)
library(extraDistr)

######################################################
########Fixed Parameters and re-used functions#######
#####################################################


#day numbering, could feed every 3 days sample mosquitoes every 7 days, sample humans every month (except sick visits)
#feed_days<- seq(from=0, to=722, by =3)
mos_sample_days<- seq(from=0,to=357,by=7)+365
human_sample_days<- c(seq(from=0,to=357,by=30),357)+365

#people
n_p<- 200



#possible  haplotypes & frequencies
haps<- rep(1:220)
freq<- c(rep(0.001,140),rep(0.002,15),rep(0.003,8),rep(0.005,4), rep(0.006,3), rep(0.007,5),
         rep(0.008,2), rep(0.009,3), rep(0.011,3), 0.012, 0.013, rep(0.015,2),rep(0.016,2),
         0.023,0.025,0.026,0.027,0.04,0.042,rep(0.044,2),0.047,0.049,0.052,0.055,0.063,0.064,0.07,
         0.077,0.084,0.093,rep(0.098,2),0.111,0.143,0.148,0.213,0.235,0.315,0.318,0.331,0.46,
         0.476,0.547)


mos_haps<- c(1:80,139:148,149,151,153,155, 156:162, 164,165,166,168,169,171,172,174,176,178,179,181,183,184,186,187,188,189,190,
             192,194,195,196,197,198,200,202,203,204,206,208,210,211,212,213,214,215,216,217,218,220) 

mos_freq<- freq[mos_haps]

human_haps<- c(60:138,143,148:155, 163:165,167,169,170,171,173,174,175,176,177,180,182,184,185,187,189,190,191,193,196,
               197,199,200,201,202,203,205,206,207,208,209,211,213,214,215,216,217,219,218,220)

human_freq<- freq[human_haps]

get_infection<- function(x){
  inf<-ifelse(x<3,rbinom(1,size=1,p=0.05),ifelse(x<4,rbinom(1,size=1,p=0.1),ifelse(x<5,rbinom(1,size=1,0.15),rbinom(1,size=1,p=0.25))))
  
  return(inf)
}
get_moi<- function(x){
  moi<- ifelse(x<3,rtpois(1,4,a=0,b=5),ifelse(x<4, rtpois(1,6,a=0,b=10), ifelse(x<5,rtpois(1,6,a=0,b=15),rtpois(1,6,a=0,b=17))))
  return(moi)
}

get_pers_infec<- function(x){
  haps_index<-sample(human_haps, size=x, prob=human_freq)
  return(haps_index)
}

get_mos_infec<- function(x){
  haps_index<-sample(mos_haps, size=x, prob=mos_freq)
  return(haps_index)
}

#haplotype exchange between human and mosquitoes if biting is happening

get_old_p_haps<- function(x){
  old_index<- which(x>=14)
  return(old_index)
}

#old_haps_p<-apply(age_haps_p,1,get_old_p_haps)

get_old_m_haps<- function(x){
  old_index<- which(x>=9)
  return(old_index)
}



#function to get the mosquito death as a function of age
get_mos_death<- function(x){
  death_m<- ifelse(x<3,rbinom(1,size=1,p=0.01),ifelse(x<6,rbinom(1,size=1,p=0.05),ifelse(x<9,rbinom(1,size=1,p=0.1),ifelse(x<12,rbinom(1,size = 1,p=0.25),ifelse(x<15,rbinom(1,size=1,p=0.5),rbinom(1,size=1,p=0.7))))))
  return(death_m)
}



#############################
####Simulation function#####
############################

run_biting_sim<- function(pr_symp, pr_clear, pr_off_feed, pr_on_feed, pr_hum_to_mos, pr_mos_to_hum,
                          pr_num_biting, n_m, scenario_name, n_sim){
  
  symptomatic_MOI_df<- matrix(NA, nrow=n_sim, ncol=n_p*722)
  mosquito_MOI_df<- matrix(NA, nrow=n_sim, ncol=30*722)
  human_MOI_df<- matrix(NA,nrow=n_sim,ncol=n_p*722)
  eir_df<- matrix(NA, nrow=n_sim,ncol=n_p)
  age_mos_df<- matrix(NA, nrow=n_sim, ncol=n_m)
  
  
  for(q in 1:n_sim){
    #mosquito ages
    
    age_m<- rtpois(n_m, 4,a=0,b=14)
    
    
    
    #starting infection status and MOIs for people
    inf_p<-rbinom(n_p, size=1, p=0.3)
    moi_p<- rep(NA, n_p)
    for(i in 1:length(inf_p)) moi_p[i]<-ifelse(inf_p[i]==1,rtpois(1,2,a=0,b=16),0)
    
    
    #starting infection status for mosquito these are mosquito-age dependent older =more likely to be infected
    
    inf_m<- sapply(age_m, get_infection)
    
    
    #starting MOIs for mosquitoes, these are also age dependent,older=higher MOI
    
    moi_m<- ifelse(inf_m==0,0,NA)
    
    
    
    moi_m[which(is.na(moi_m))]<-sapply(age_m[which(is.na(moi_m))],get_moi)
    
    #starting haplotypes for people and age of those haplotypes
    
    
    infec_p<- sapply(moi_p,get_pers_infec)
    
    
    
    age_haps_p<- matrix(0,nrow=n_p, ncol=length(haps))
    
    for(i in 1:n_p){
      if(length(infec_p[[i]])>0){
        age_haps_p[i,infec_p[[i]]]<- rpois(length(infec_p[[i]]),20)
      }
      
    }
    #starting haplotypes for mosquitoes  and age of those haplotypes
    
    
    
    
    infec_m<- sapply(moi_m,get_mos_infec)
    
    
    age_haps_m<- matrix(0,nrow=n_m, ncol=length(haps))
    
    for(i in 1:n_m){
      if(length(infec_m[[i]])>0){
        age_haps_m[i,infec_m[[i]]]<- sample(1:age_m[i],size=length(infec_m[[i]]), replace=T)
      }
      
    }
    
    
    
    
    symp_index<- rep(0,n_p)
    
    symp_index[sample(which(inf_p==1),size=floor(length(which(inf_p==1))*pr_symp))]<-1
    
    symp_age<- rep(0,n_p)
    symp_age[which(symp_index==1)]<- 1
    
    
    num_infec_bites<- rep(0,n_p)
    
    
    
    ##############################
    ######Start of timed sim#####
    
    for(r in 1:722){
      
      #age 1 day 
      age_m<- age_m+1
      
      age_haps_m[age_haps_m!=0]<- age_haps_m[age_haps_m!=0]+1
      
      age_haps_p[age_haps_p!=0]<- age_haps_p[age_haps_p!=0]+1
      #is the mosquito dead
      death_m<- sapply(age_m,get_mos_death)
      #if so replace it with a new mosquito(newborn)
      
      age_m[which(death_m==1)]<- 0
      
      
      inf_m[which(death_m==1)]<- 0
      
      
      moi_m[which(death_m==1)]<- 0
      
      infec_m[which(death_m==1)]<- NA
      
      for(i in 1:n_m){
        if(death_m[i]==1 ){
          
          
          age_haps_m[i,]<- 0
          
          
          
        }
      }
      
      
      
      # symptomatic infections get  sampled (sick visit) and cleared immediately (drugs)
      
      symp_age[which(symp_age!=0)]<- symp_age[which(symp_age!=0)]+1
      symp_index<- ifelse(symp_age>14,0,symp_index)
      symp_age<- ifelse(symp_age>14,0,symp_age)
      
      infec_p[which(symp_age==2)]<- NA
      age_haps_p[which(symp_age==2),]<-0
      
      
      for(i in 1:n_p){
        if(length(which(age_haps_p[i,]>=7))>0& length(which(age_haps_p[i,]<=30))>0&max(age_haps_p[i,])!=0
           &symp_index[i]==0){
          new_index<-rbinom(1,1,pr_symp)
          if(new_index==1){
            symp_index[i]<-1
            symp_age[i]<- 1
          }
        }
      }
      
      symptomatics<- which(symp_age==1)
      if(length(symptomatics)>0){
        symp_moi<- rep(NA, length(symptomatics))
        for(i in 1:length(symptomatics)){
          symp_moi[i]<- length(na.omit(unlist(infec_p[[symptomatics[i]]])))
        }
        if(sum(symp_moi,na.rm=T)>0){
          symptomatic_MOI_df[q,(1+(n_p*(r-1))):(((n_p*(r-1)))+length(symp_moi))]<- symp_moi
        }
      }
      
      
      
      
      
      
      
      #then clear human parasites by age of parasites
      
      
      
      for(i in 1:n_p){
        if(length(na.omit(infec_p[[i]]))>0){
          for(j in c(na.omit(infec_p[[i]]))){
            
            clear<- ifelse(age_haps_p[i,j]>=90,rbinom(1,1,0.95),ifelse(age_haps_p[i,j]>=65,rbinom(1,1,0.85),ifelse(age_haps_p[i,j]>=30,rbinom(1,1,pr_clear),0)))
          }
          
          if(clear==1){
            infec_p[[i]][which(infec_p[[i]]==j)]<- NA
            age_haps_p[i,j]<- 0
          }
        }
      }
      
      
      
      
      
      
      #if its a feeding day
      
      
      
      mos_bite<- matrix(0, n_m,n_p)
      for(i in 1:n_m){
        if(age_m[i]>=2){
          
          bite<- ifelse(min(age_haps_m[i,])<=(age_m[i]-3),rbinom(1,1,pr_off_feed),rbinom(1,1,pr_on_feed))
          
          if(bite==1){
            num_biting<- sample(c(1,2,3,4,5,6,7),size=1,prob=pr_num_biting)
            
            mos_bite[i,sample(1:200,size=num_biting,replace=F)]<-1
          }
        }
        
      }
      
      #haplotype exchange if biting is happening
      
      old_haps_p<-apply(age_haps_p,1,get_old_p_haps)
      
      
      old_haps_m<-apply(age_haps_m,1,get_old_m_haps)
      
      ##########################
      which_mos_bite<- ifelse(apply(mos_bite,1,sum)>0,1,0)
      
      old_hap_mosquitoes<- rep(0,n_m)
      
      for(i in 1:n_m){
        if(length(old_haps_m)>0){
          if(length(old_haps_m[[i]])>0){
            old_hap_mosquitoes[i]<-1
          }
        }
      }
      
      hap_transfer_mos<- rep(0,n_m)
      
      for(i in 1:n_m){
        if(old_hap_mosquitoes[i]==1&which_mos_bite[i]==1){
          hap_transfer_mos[i]<-1
        }
      }
      
      
      ####################
      
      which_hum_bite<- ifelse(apply(mos_bite,2,sum)>0,1,0)
      
      old_hap_hum<- rep(0,n_p)
      
      
      for(i in 1:n_p){
        if(length(old_haps_p)>0){
          if(length(old_haps_p[[i]])>0){
            old_hap_hum[i]<-1
          }
        }
      }
      
      hap_transfer_hum<- rep(0,n_p)
      
      for(i in 1:n_p){
        if(old_hap_hum[i]==1&which_hum_bite[i]==1){
          hap_transfer_hum[i]<-1
        }
      }
      
      ##########################
      
      
      
      for(i in 1:n_p){
        mos_index<- which(mos_bite[,i]==1)
        if(length(mos_index)>0){
          inf_bites<- rep(0,length(mos_index))
          for(j in 1:length(mos_index)){
            if(length(old_haps_p)>0){
              if(length(old_haps_p[[i]])>0&symp_age[i]<2){
                transfer_haps<- rep(NA,length(old_haps_p[[i]]))
                for(k in 1:length(old_haps_p[[i]])){
                  prob<- pr_hum_to_mos
                  transfer<- rbinom(1,1,prob)
                  if(transfer==1){
                    transfer_haps[k]<- old_haps_p[[i]][k]
                  }
                }
                new_haps<- setdiff(na.omit(transfer_haps),infec_m[[mos_index[j]]])
                infec_m[[mos_index[j]]]<-c(infec_m[[mos_index[j]]],new_haps)
                age_haps_m[mos_index[j],new_haps]<-1
              }
            }
            
            if(length(old_haps_m)>0){
              if(length(old_haps_m[[mos_index[j]]])>0&symp_index[i]<1){
                transfer_haps_m<- rep(NA,length(old_haps_m[[mos_index[j]]]))
                for(l in 1:length(old_haps_m[[mos_index[j]]])){
                  prob_m<- pr_mos_to_hum
                  transfer_m<- rbinom(1,1,prob_m)
                  if(transfer_m==1){
                    transfer_haps_m[l]<- old_haps_m[[mos_index[j]]][l]
                  }
                }
                new_haps_m<- setdiff(na.omit(transfer_haps_m), infec_p[[i]])
                infec_p[[i]]<- c(infec_p[[i]],new_haps_m)
                age_haps_p[i,new_haps_m]<- 1
                inf_bites[j]<-1
              }
              
              
              
            }
          }
          
          num_infec_bites[i]<- num_infec_bites[i]+sum(inf_bites)
        }
        
      }
      
      
      #if its a human and mosquito sample day     
      if(r%in%human_sample_days&r%in%mos_sample_days){
        sample_index<- sample(1:n_m,size=30)
        mos_moi<- rep(NA,30)
        for(t in 1:30){
          mos_moi[t]<- length(na.omit(unlist(infec_m[[sample_index[t]]])))
        }
        mosquito_MOI_df[q,(1+(30*(r-1))):(30*r)]<- mos_moi
        
        
        
        age_m[sample_index]<- 0
        #rtpois(length(sample_index), 4,a=1,b=14)
        
        inf_m[sample_index]<- 0
        #sapply(age_m[sample_index],get_infection)
        
        moi_m[c(sample_index)]<- 0
        #,which(inf_m==1)
        #sapply(age_m[c(sample_index,which(inf_m==1))],get_moi)
        infec_m[c(sample_index)]<-0   
        #,which(inf_m==1&moi_m>0)
        #sapply(moi_m[c(sample_index,which(inf_m==1&moi_m>0))],get_mos_infec)
        for(i in 1:n_m){
          if(i%in%sample_index ){
            #& inf_m[i]==1 & moi_m[i]>0
            
            age_haps_m[i,]<- 0
            #c(infec_m[[i]])
            #sample(1:age_m[i],size=length(infec_m[[i]]),replace=T)
            
          }
        }
        hum_moi<- rep(NA, n_p)
        for(w in 1:n_p){
          hum_moi[w]<- length(na.omit(unlist(infec_p[[w]])))
        }
        
        human_MOI_df[q,((1+(200*(r-1))):(200*r))]<- hum_moi
        
      }  
      #if its only a mosquito sample day
      if(r%in% mos_sample_days){
        sample_index<- sample(1:n_m,size=30)
        mos_moi<- rep(NA,30)
        for(t in 1:30){
          mos_moi[t]<- length(na.omit(unlist(infec_m[[sample_index[t]]])))
        }
        mosquito_MOI_df[q,(1+(30*(r-1))):(30*r)]<- mos_moi
        
        
        
        age_m[sample_index]<- 0
        #rtpois(length(sample_index), 4,a=1,b=14)
        
        inf_m[sample_index]<- 0
        #sapply(age_m[sample_index],get_infection)
        
        moi_m[c(sample_index)]<- 0
        #,which(inf_m==1)
        #sapply(age_m[c(sample_index,which(inf_m==1))],get_moi)
        infec_m[c(sample_index)]<- NA
        #,which(inf_m==1&moi_m>0)
        #sapply(moi_m[c(sample_index,which(inf_m==1&moi_m>0))],get_mos_infec)
        for(i in 1:n_m){
          if(i%in%sample_index ){
            
            #& inf_m[i]==1 & moi_m[i]>0
            age_haps_m[i,]<- 0
            #c(infec_m[[i]])
            #sample(1:age_m[i],size=length(infec_m[[i]]),replace=T)
          }
        }
        
        
      }
      #if its only a human sample day
      if(r%in% human_sample_days){
        
        hum_moi<- rep(NA, n_p)
        for(w in 1:n_p){
          hum_moi[w]<- length(na.omit(unlist(infec_p[[w]])))
        }
        
        human_MOI_df[q,((1+(200*(r-1))):(200*r))]<- hum_moi
        
        
      }
      
    }
    
    eir_df[q,]<- num_infec_bites
    age_mos_df[q,]<- age_m
    print(q)
  }
  saveRDS(symptomatic_MOI_df,file=paste0("symptomatic_MOI_",scenario_name))
  saveRDS(mosquito_MOI_df, file =paste0("mosquito_MOI_",scenario_name))
  saveRDS(human_MOI_df,file =paste0("human_MOI_",scenario_name))
  saveRDS(eir_df, file =paste0("eir_", scenario_name))
  saveRDS(age_mos_df, file=paste0("mos_age_",scenario_name))
}

run_biting_sim(pr_symp = 0.1, pr_clear = 0.75, pr_off_feed = 0.04, pr_on_feed = 0.4,pr_hum_to_mos = 1.0, pr_mos_to_hum = 0.65, pr_num_biting = c(0.025,0.025,0.1,0.7,0.1,0.025,0.025), n_m=30000, scenario_name = "BitingA", n_sim=50)
run_biting_sim(pr_symp = 0.1, pr_clear = 0.75, pr_off_feed = 0.04, pr_on_feed = 0.4,pr_hum_to_mos = 1.0, pr_mos_to_hum = 0.65, pr_num_biting = c(0.7,0.1,0.1,0.025,0.025,0.025,0.025), n_m=30000, scenario_name = "BitingB", n_sim=50)
run_biting_sim(pr_symp = 0.1, pr_clear = 0.75, pr_off_feed = 0.04, pr_on_feed = 0.4,pr_hum_to_mos = 1.0, pr_mos_to_hum = 0.65, pr_num_biting = c(0.6,0.35,0.04,0.01,0,0,0), n_m=30000, scenario_name = "BitingC", n_sim=50)
run_biting_sim(pr_symp = 0.1, pr_clear = 0.75, pr_off_feed = 0.04, pr_on_feed = 0.4,pr_hum_to_mos = 1.0, pr_mos_to_hum = 0.65, pr_num_biting = c(0.65,0.35,0,0,0,0,0), n_m=30000, scenario_name = "BitingD", n_sim=50)
run_biting_sim(pr_symp = 0.1, pr_clear = 0.75, pr_off_feed = 0.04, pr_on_feed = 0.4,pr_hum_to_mos = 1.0, pr_mos_to_hum = 0.65, pr_num_biting = c(1,0,0,0,0,0,0), n_m=30000, scenario_name = "BitingE", n_sim=50)



