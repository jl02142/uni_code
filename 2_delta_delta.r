# J Ludwig 2021

# Delta Delta method for real time PCR to normalize and calculate fold change

    getwd()
    setwd()
    input <- read.csv("input2.csv", row.names=1)
    

    Tot_colm=length(input)
    Tot_samp=length(input[,1])
    Half_samp=Tot_samp/2

    Avg_cont_HKG=mean(input[,1][1:Half_samp])    
    Avg_expi_HKG=mean(input[,1][(Half_samp+1):Tot_samp])

    Avg_cont_GOI <- array(-9999, dim=c(1, Tot_colm))
    Avg_expi_GOI <- array(-9999, dim=c(1, Tot_colm))
    Delta_CT_cont <- array(-9999, dim=c(1, Tot_colm))
    Delta_CT_expi <- array(-9999, dim=c(1, Tot_colm))
    Delta_delta_CT <- array(-9999, dim=c(1, Tot_colm))
    Expr_fold_chng <- array(-9999, dim=c(1, Tot_colm))

    for(i in 2:Tot_colm){
        Avg_cont_GOI[i]=mean(input[,i][1:Half_samp])
        Avg_expi_GOI[i]=mean(input[,i][(Half_samp+1):Tot_samp])        
    }

    for(i in 2:Tot_colm){
        Delta_CT_cont[i] = Avg_cont_GOI[i] - Avg_cont_HKG
        Delta_CT_expi[i] = Avg_expi_GOI[i] - Avg_expi_HKG
    }

    for(i in 2:Tot_colm){
        Delta_delta_CT[i] = Delta_CT_expi[i] - Delta_CT_cont[i]
        Expr_fold_chng[i] = 2^-Delta_delta_CT[i]
    }

    #Avg_cont_GOI=mean(input[,2][1:Half_samp])
    #Avg_expi_GOI=mean(input[,2][(Half_samp+1):Tot_samp])

    #Delta_CT_cont=Avg_cont_GOI - Avg_cont_HKG
    #Delta_CT_expi=Avg_expi_GOI - Avg_expi_HKG

    #Delta_delta_CT=Delta_CT_expi - Delta_CT_cont
    #Expr_fold_chng=2^-Delta_delta_CT

    Indi_fold_change <- array(-9999, dim=c(Tot_samp, Tot_colm))

    #Indi_fold_change<- vector(mode = "numeric", Tot_samp)
    for(i in 2:Tot_colm){
        for(j in 1:Tot_samp){
            Indi_fold_change[,i][j] = 2^-((input[,i][j] - input[,1][j]) - Delta_CT_cont[i])
        }
    }
    
    Meth_1_cont_avg <- array(-9999, dim=c(1, Tot_colm))
    Meth_1_expi_avg <- array(-9999, dim=c(1, Tot_colm))
    Meth_1_cont_std <- array(-9999, dim=c(1, Tot_colm))
    Meth_1_expi_std <- array(-9999, dim=c(1, Tot_colm))

    for(i in 2:Tot_colm){

            Meth_1_cont_avg[i] = mean(Indi_fold_change[,i][1:Half_samp])
            Meth_1_expi_avg[i] = mean(Indi_fold_change[,i][(Half_samp+1):Tot_samp])
            Meth_1_cont_std[i] = sd(Indi_fold_change[,i][1:Half_samp])
            Meth_1_expi_std[i] = sd(Indi_fold_change[,i][(Half_samp+1):Tot_samp])

    }
    
    #Meth_1_cont_avg = mean(Indi_fold_change[1:Half_samp])
    #Meth_1_expi_avg = mean(Indi_fold_change[(Half_samp+1):Tot_samp])
    #Meth_1_cont_std = sd(Indi_fold_change[1:Half_samp])
    #Meth_1_expi_std = sd(Indi_fold_change[(Half_samp+1):Tot_samp])

    Meth_1_cont_SE <- array(-9999, dim=c(1, Tot_colm))
    Meth_1_expi_SE <- array(-9999, dim=c(1, Tot_colm))

    #Meth_1_cont_SE = Meth_1_cont_std / sqrt(Half_samp)
    #Meth_1_expi_SE = Meth_1_expi_std / sqrt(Half_samp)

    Data_for_ttest <- array(-9999, dim=c(Tot_samp, Tot_colm))
    for(i in 2:Tot_colm){
        Meth_1_cont_SE[i] = Meth_1_cont_std[i] / sqrt(Half_samp)
        Meth_1_expi_SE[i] = Meth_1_expi_std[i] / sqrt(Half_samp)
        for(j in 1:Tot_samp){
            Data_for_ttest[,i][j] = input[,i][j]-input[,1][j]
        }
    }

    stor_ttest <- array(-9999, dim=c(1, Tot_colm))

    for(i in 2:Tot_colm){
        print(t.test(Data_for_ttest[,i][1:Half_samp], Data_for_ttest[,i][(Half_samp+1):Tot_samp], var.equal=TRUE))
        stor_ttest[i] = t.test(Data_for_ttest[,i][1:Half_samp], Data_for_ttest[,i][(Half_samp+1):Tot_samp], var.equal=TRUE)$p.value
    }
    
    col_names_print=colnames(input[2:Tot_colm])
    df_print<-data.frame(matrix(nrow=Tot_colm-1,ncol=5, dimnames=list(c(col_names_print), c("p-value", "cont_mean", "cont_se", "expi_mean", "expi_se"))))

    for(i in 2:Tot_colm){
        df_print[,1][i-1] = stor_ttest[i]
        df_print[,2][i-1] = Meth_1_cont_avg[i]
        df_print[,3][i-1] = Meth_1_cont_SE[i]
        df_print[,4][i-1] = Meth_1_expi_avg[i]
        df_print[,5][i-1] = Meth_1_expi_SE[i]
    }

    df_print
    
    write.csv(df_print, file=paste0("2_delta_delta_stats_", format(Sys.time(), "%m_%d_%H_%M_%S"), ".csv"))
