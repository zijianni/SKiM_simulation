library(qvalue)

get_OR <- function(p_1, p_2){
    # Calculate odds ratio between p_1 and p_2
    # OR=[p_1/(1-p_1)]/[p_2/(1-p_2)]
    # Args:
    #   p_1 (num): P(B occurs | A occurs)
    #   p_2 (num): P(B occurs | A not occurs)
    # Returns: 
    #   (num) odds ratio
    
    OR <- (p_1/(1-p_1))/(p_2/(1-p_2))
    return(OR)
}

get_p1 <- function(p_2, OR){
    # Given p_2 and odds ratio, calculate p_1
    # Args:
    #   p_2 (num): P(B occurs | A not occurs)
    #   OR (num): odds ratio between p_1 and p_2
    # Returns: 
    #   (num) P(B occurs | A occurs)
    
    p_1 <- 1/(1+(1-p_2)/p_2/OR)
    return(p_1)
}


simulate_table <- function(n_paper,
                           n_paper_have_x, 
                           n_paper_no_x=NULL,
                           p_y_in_paper_have_x, 
                           p_y_in_paper_no_x){
    # Simulate contingency table
    # Args:
    #   n_paper (int): total number of papers
    #   n_paper_have_x (int): number of papers containing term x
    #   n_paper_no_x (int): number of papers not containing term x
    #   p_y_in_paper_have_x (num): 
    #       probability of observing term y in papers containing term x
    #   p_y_in_paper_no_x (num): 
    #       probability of observing term y in papers not containing term x
    # Returns: 
    #   (matrix) a 2x2 contingency table
    
    if(is.null(n_paper_no_x)){
        n_paper_no_x <- n_paper-n_paper_have_x
    }
    if((n_paper_have_x+n_paper_no_x!=n_paper)){
        stop("Invalid input values.")
    }
    
    x1 = rbinom(1, n_paper_have_x, p_y_in_paper_have_x)
    x2 = rbinom(1, n_paper_no_x, p_y_in_paper_no_x)
    have_x = c(x1, n_paper_have_x - x1)
    no_x = c(x2, n_paper_no_x - x2)
    TBL = rbind(have_x, no_x)#+1 # add pseudo-counts to prevent 0's
    colnames(TBL) <- c("have_y","no_y")
    
    return(TBL)
}

p_adjust <- function(pvalue){
    # Perform FDR control of multiple testing using BH, q-value and local FDR
    # Args:
    #   pvalue (vector of num): p-values from multiple testing
    # Returns:
    #   (data.frame) columns: 
    #       p_val (num): original input p-values
    #       p_adj_<method> (num): adjusted p-values under each method
    
    pvalue <- pmin(pvalue, 1)
    
    # Regular FDR control
    df_out <- data.frame(p_val=pvalue)
    
    df_out$p_adj_BH <- p.adjust(pvalue, method = "BH")
    
    qobj <- qvalue(p = pvalue)
    df_out$p_adj_qvalue <- qobj$qvalues
    df_out$p_adj_locfdr <- qobj$lfdr
    
    return(df_out)
    
}

P_adjust <- function(pvalue_mat){
    # Perform FDR control of matrix of p-values.
    # Designed for all B->C pairs
    # Args:
    #   pvalue_mat (matrix of num): p-values. Rows are Cs. Columns are Bs.
    # Returns:
    #   (list):
    #       p_val (matrix of num): original input p-values
    #       p_adj_<method> (matrix of num): adjusted p-values under each method
    
    n_row <- nrow(pvalue_mat)
    n_col <- ncol(pvalue_mat)

    p_adj <- p_adjust(as.vector(pvalue_mat))
    
    P_adj <- list()
    for(nm in colnames(p_adj)){
        P_adj[[nm]] <- matrix(p_adj[,nm],n_row, n_col)
    }

    return(P_adj)
}

get_power_FDR <- function(predicted_label, true_label){
    # Calculate power and FDR of certain prediction
    # Args:
    #   predicted_label (vector of binary or logical): predicted significancy. 
    #       1 (TRUE): significant. 0 (FALSE): not significant.
    #   true_label (vector of binary or logical): true significancy. 
    #       1 (TRUE): significant. 0 (FALSE): not significant.
    # Returns:
    #   (list):
    #       power: proportion of predicted significant events
    #           among all true significant events.
    #       FDR: proportion of true insignificant events
    #           among all predicted significant events
    predicted_label <- as.integer(predicted_label)
    true_label <- as.integer(true_label)
    
    true_significant <- true_label==1
    predicted_significant <- predicted_label==1
    power <- mean(predicted_label[true_significant])
    FDR <- 1-mean(true_label[predicted_significant])
    
    return(list(power=power, FDR=FDR))
}

test_CtoB <- function(cand_B, p_C, which_signif_CtoB, verbose=TRUE){
    # Simulate contingency table and perform FET between B's and C's
    # Args:
    #   cand_B (vector of int): indices of candidate B's to test
    #   p_C (vector of num): estimated probability of occurrence for each C
    #   which_signif_CtoB (matrix of int): the significancy matrix between Bs and Cs
    #   verbose (logical): whether or not print current progress.
    # Returns:
    #   (list):
    #       PVAL_CtoB (matrix of num): p-value matrix
    #       SORT_RATIO_CtoB (matrix of num): sort ratio matrix
    #       SCORE_CtoB (matrix of num): prediction score matrix
    
    PVAL_CtoB <- matrix(0,length(p_C),length(cand_B))
    colnames(PVAL_CtoB) <- cand_B
    SORT_RATIO_CtoB <- SCORE_CtoB <- PVAL_CtoB
    
    for(i in seq_along(cand_B)){
        B_idx <- cand_B[i]
        if(verbose){
            message("Currently calculating B_",B_idx,".")
        }
        
        # Set p_1 and p_2 for B and each C
        p_1_CtoB <- p_2_CtoB <- p_C
        OR_CtoB <- 2+rexp(nrow(which_signif_CtoB))
        
        # Get the 10 significant C's for this B
        which_signif_C <- which_signif_CtoB[,B_idx]
        
        # # We choose a large p_2 since FET p-values are discrete and couldn't distinguish too tiny p_1 and p_2
        # p_1_CtoB[which_signif_C] <- pmax(get_p1(p_2_CtoB[which_signif_C],OR_CtoB),1e-2)
        
        # Simulate contingency tables and perform FET between B and each C
        pval_CtoB <- sort_ratio_CtoB <- numeric(nrow(meta_C))
        
        for(C_idx in seq_along(pval_CtoB)){
            temp_table <- simulate_table(N, n_B[B_idx], N-n_B[B_idx], p_1_CtoB[C_idx], p_2_CtoB[C_idx])
            
            # Add pseudocounts for associated B->C pairs
            if(C_idx%in%which_signif_C){
                pseudo_mat <- min(5,sum(temp_table[,1]),sum(temp_table[1,]))*matrix(c(1,-1,-1,1),2,2)
                temp_table <- temp_table+pseudo_mat
            }
            
            pval_CtoB[C_idx] <- fisher.test(temp_table)$p.value
            sort_ratio_CtoB[C_idx] <- temp_table[1,1]/sum(temp_table[,1])
        }
        PVAL_CtoB[,i] <- pval_CtoB
        SORT_RATIO_CtoB[,i] <- sort_ratio_CtoB
        SCORE_CtoB[,i] <- -log10(pval_CtoB)+log10(sort_ratio_CtoB)
        # For groups with small n_B[B_idx], p-values is not meaningful as well.
        # I will expect this happen in practice in SKiM. i.e.
        # SKiM's findings are more focused on categories with more occurrence.
        
        # # Calculate adjusted p-values
        # padj_CtoB <- p_adjust(pval_CtoB)
        
    }
    return(list(PVAL_CtoB=PVAL_CtoB,
                SORT_RATIO_CtoB=SORT_RATIO_CtoB,
                SCORE_CtoB=SCORE_CtoB))
}

