###################
# Zijian Ni
# 10/01/2020

set.seed(2020)  # for reproducibility
SEEDS <- rnorm(100)
setwd("~/Google Drive/Hallu/codes/ckgroup/SKIM/")
source("simulation_functions.R")

N <- 29613663 # total number of papers in database
N_A <- c(500, 5000, 50000, 100000, 200000) # number of papers containing A
meta_B <- readxl::read_xlsx("SKiM_Files/To_Zijian/Phenotypes_and_symptoms_count.xlsx", 
                            sheet = 1)
meta_C <- readxl::read_xlsx("SKiM_Files/To_Zijian/Drugs_count.xlsx",
                            sheet = 1)
n_B <- meta_B$Phenotype_and_symptom_count # number of papers containing each B
n_C <- meta_C$Drug_count # number of papers containing each C

# Filter out zero occurrence terms
meta_B <- meta_B[n_B>0,]
n_B <- n_B[n_B>0]
meta_C <- meta_C[n_C>0,]
n_C <- n_C[n_C>0]

p_B <- n_B/N
p_C <- n_C/N

n_signif_B <- 20 # number of significant B's
n_signif_C <- 10 # number of significant C's for each B

# For significant B's in A->B FET, keep top_N B's 
# with highest prediction score for B->C tests.
top_N <- 50 

n_rep <- 10 # number of replicates under each setting




##########################
# For each A, repeat n_rep times 
#########################

Res <- list() # Results of all A's

for(a_idx in seq_along(N_A)){
    message("Now working on ",a_idx,"/",length(N_A)," of A's")
    n_A <- N_A[a_idx]
    Res_tmp <- list() # Results of current A
    set.seed(SEEDS[a_idx])
    for(rp in seq_len(n_rep)){
        message("Now working on ",rp,"/",n_rep, " of replicates.")

        #######################################
        # Set significant B's
        which_signif_BtoA <- sample.int(nrow(meta_B), n_signif_B) 
        
        #######################################
        # Set significant C's for each B
        which_signif_CtoB <- replicate(nrow(meta_B),sample.int(nrow(meta_C), n_signif_C))
        
        #######################################
        # Set significant C's to A
        which_signif_CtoA <- unique(as.vector(which_signif_CtoB[,which_signif_BtoA]))
        
        
        #######################################
        # Set p_1 and p_2 for A and each B
        p_1_BtoA <- p_2_BtoA <- p_B
        OR_BtoA <- 2+rexp(n_signif_B)
        p_1_BtoA[which_signif_BtoA] <- get_p1(p_2_BtoA[which_signif_BtoA],OR_BtoA)
        
        
        #######################################
        # Simulate contingency tables, calculate sort ratios and perform FET between A and each B
        pval_BtoA <- sort_ratio_BtoA <- numeric(nrow(meta_B))
        
        for(B_idx in seq_along(pval_BtoA)){
            temp_table <- simulate_table(N, n_A, N-n_A, 
                                         p_1_BtoA[B_idx], p_2_BtoA[B_idx])
            # Add pseudo-count for true associations
            if(B_idx%in%which_signif_BtoA){
                pseudo_mat <- min(5,temp_table[1,2],temp_table[2,1])*
                    matrix(c(1,-1,-1,1),2,2)
                temp_table <- temp_table+pseudo_mat
            }
            sort_ratio_BtoA[B_idx] <- temp_table[1,1]/sum(temp_table[,1])
            pval_BtoA[B_idx] <- fisher.test(temp_table)$p.value
        }
        
        # Calculate prediction score
        score_BtoA <- -log10(pval_BtoA)+log10(sort_ratio_BtoA)
        
        ######################################
        # Keep significant Bs with p-value <=1e-5, then
        # choose top_N B's with largest prediction score
        # When there are less than top_N significant B's,
        # just use all these significant B's.
        signif_B <- which(pval_BtoA<=1e-5)
        cand_B <- signif_B[rank(-score_BtoA[signif_B])<=top_N]

        
        ####################################
        # Simulate contingency table and perform FET for top $N$ significant $B$'s and all $C$'s
        res_CtoB <- test_CtoB(cand_B,p_C,which_signif_CtoB, verbose=T)
        
        
        ####################################
        # Calculate power and FDR for BtoA
        p_F_BtoA <- get_power_FDR(pval_BtoA<=1e-5, 
                                  seq_along(pval_BtoA)%in%which_signif_BtoA)

        ####################################
        # Calculate power and FDR for CtoB
        # True significancy
        true_signif_CtoB <- apply(which_signif_CtoB[,cand_B],
                                  2, function(x) seq_len(nrow(meta_C))%in%x)
        p_F_CtoB <- get_power_FDR(as.vector(res_CtoB$PVAL_CtoB<=1e-5), 
                                  as.vector(true_signif_CtoB))
        
        ####################################
        # Calculate power and FDR for CtoA
        SKiM_signif_CtoA <- apply(res_CtoB$PVAL_CtoB,1,function(x) any(x<=1e-5))
        
        p_F_CtoA <- get_power_FDR(SKiM_signif_CtoA, 
                                  seq_len(nrow(meta_C))%in%which_signif_CtoA)
        
        
        ###################################
        # Write results
        Res_tmp[[rp]] <- list(res_CtoB=res_CtoB, 
                              res_BtoA=list(
                                  pval_BtoA=pval_BtoA, 
                                  sort_ratio_BtoA=sort_ratio_BtoA, 
                                  score_BtoA=score_BtoA
                                  ), 
                              p_F_CtoB=p_F_CtoB, 
                              p_F_BtoA=p_F_BtoA, 
                              p_F_CtoA=p_F_CtoA)
    }
    
    names(Res_tmp) <- paste0("Rep_",seq_len(n_rep))
    Res[[a_idx]] <- Res_tmp
}

names(Res) <- paste0("n_A=",N_A)

Res$session <- sessionInfo()

save(Res, file=paste0("sim_out_",Sys.Date(),".RData"))


