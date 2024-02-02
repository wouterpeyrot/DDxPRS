show("Load libraries: mvtnorm, mvnfast")

DDxPRS <- function( prs_liab , prs_r2l , snp_h2l , K  , crosstrait_cor.prs , crosstrait_rg, clinical.prior , liab.configuration  ){

	## Before applying DDxPRS, it is highly encouraged to study and run the example at https://github.com/wouterpeyrot/DDxPRS 

	## All disorder names should be provided in every input item. Disorders should be identically ordered for all input items.
	## prs_liab: dataframe with the liability-scale prs for every disorder
	## prs_r2l: vector with liability-scale variance explained in every disorder by its case-control prs
	## snp_h2l: vector with liability-scale heritabilities for every disorder
	## K: vector with population prevalence for every disorder
	## crosstrait_cor.prs: dataframe with correlation between the case-control prs in the population reference sample
	## crosstrait_rg: datafrane with genetic correlation estimated between the disorders
	## clinical.prior: a vector with prior clinical probability the diagnostic categories defined in 
	## liab.configuration: - A dataframe linking the configurations of liabilities with the diagnostic categories. 
	##                     - Names of the first columns of liab.configuration should match the disorders
	##                     - The last columns should be named "diagnostic.category"

	#####################################
	## double-check input
	#####################################

		disorder_names <- names(K)
		n_disorders <- length(K)
		category_names <- names(clinical.prior)
		analytical_col <- disorder_names
		liab.configuration$quick_name <- apply( liab.configuration[,analytical_col], MAR=1, FUN=function(x){paste(x,collapse="")} )

		if( !all.equal(names(prs_r2l),names(snp_h2l))==TRUE                ){ stop("Disorder names in prs_r2l and snp_h2l should exactly match, and identically ordered.")  } 
		if( !all.equal(names(prs_r2l),names(K))==TRUE                      ){ stop("Disorder names in prs_r2l and K should exactly match, and identically ordered.")  } 
		if( !all.equal(names(prs_r2l),colnames(crosstrait_cor.prs))==TRUE  ){ stop("Disorder names in prs_r2l and colnames(crosstrait_cor.prs) should exactly match, and identically ordered.")  } 
		if( !all.equal(names(prs_r2l),rownames(crosstrait_cor.prs))==TRUE  ){ stop("Disorder names in prs_r2l and rownames(crosstrait_cor.prs) should exactly match, and identically ordered.")  } 
		if( !all.equal(names(prs_r2l),colnames(crosstrait_rg))==TRUE       ){ stop("Disorder names in prs_r2l and colnames(crosstrait_rg) should exactly match, and identically ordered.")  } 
		if( !all.equal(names(prs_r2l),rownames(crosstrait_rg))==TRUE       ){ stop("Disorder names in prs_r2l and rownames(crosstrait_rg) should exactly match, and identically ordered.")  } 
		if( sum({names(clinical.prior) %in% liab.configuration$diagnostic.category}==FALSE) > 0 ){ stop("All diagnostic category names in clinical.prior should appear in liab.configuration$diagnostic.category.")  } 
		if( sum(clinical.prior)!=1 ){ stop("The clinical prior probabilities should add up to 1.")  } 
		if( dim(prs_liab)[2]!=n_disorders ){ stop("prs_liab should contain exactly one prs per disorder, identically ordered to prs_r2l, snp_h2l, K and the rows in crosstrait_cor.prs, crosstrait_rg")  }
		if( !all.equal(colnames( liab.configuration )[1:n_disorders],names(prs_r2l))==TRUE  ){ stop("The first columns of liab.configuration should exactly match the disorder names in prs_r2l, and identically ordered.")  } 

	#####################################
	## define several expected properties
	#####################################

		## add analyt_pop_proportion to liab.configuration
		varcov_g    <- crosstrait_rg * sqrt( cbind(snp_h2l) %*% rbind(snp_h2l) )
		varcov_env  <- crosstrait_rg * sqrt( cbind(1-snp_h2l ) %*% rbind(1-snp_h2l ) )
		varcov_liab <- varcov_g + varcov_env ## only very slight deviation from crosstrait_rg
		rownames(varcov_liab) <- colnames(varcov_liab)
		thresholds  <- -qnorm(K,mean=0,sd=1) 
		for(row in 1:dim(liab.configuration)[1]){
			temp <- as.data.frame( array(NA,dim=c(n_disorders,2)) ); colnames(temp)<-c("lower","upper") 
			for(j in 1:n_disorders ){
				temptemp <- liab.configuration[row,analytical_col[j]]
				if(temptemp==1){ temp$lower[j]<-thresholds[j] ; temp$upper[j]<-Inf   }
				if(temptemp==0){ temp$lower[j]<- -Inf         ; temp$upper[j]<-thresholds[j]   }
			}
			liab.configuration$analyt_pop_proportion[row] <- pmvnorm( lower=temp$lower , upper=temp$upper , sigma=as.matrix(varcov_liab) )[1]
		}

		## add test_priorprob to liab.configuration
		liab.configuration$test_priorprob <- NA
		for(category_name in category_names){
			index <- liab.configuration$diagnostic.category==category_name
			liab.configuration$test_priorprob[index] <- clinical.prior[category_name] * liab.configuration$analyt_pop_proportion[index] / sum(liab.configuration$analyt_pop_proportion[index])
		}

	#####################################
	## expected VARCOV of prs and liability in full population
	#####################################

		cov_prsAliabB_theory <- function( rg_AB, prs_r2l_A,h2l_A,h2l_B  ){
			sqrt( (prs_r2l_A/h2l_A) * rg_AB^2)*sqrt(prs_r2l_A*h2l_B) 
		}

		colnames <- paste( rep(c("prs_","liab_"),each=n_disorders),1:n_disorders,sep=""  )
		VARCOV <- as.data.frame(array(NA,dim=c(length(colnames),length(colnames)))) ; colnames(VARCOV)<-rownames(VARCOV)<-colnames
		for(i in 1:n_disorders){
			VARCOV[paste("prs_",i,sep="") ,paste("prs_",i,sep="") ]  <- prs_r2l[i]
			VARCOV[paste("prs_",i,sep="") ,paste("liab_",i,sep="") ] <- prs_r2l[i]
			VARCOV[paste("liab_",i,sep=""),paste("prs_",i,sep="") ]  <- prs_r2l[i]
			VARCOV[paste("liab_",i,sep=""),paste("liab_",i,sep="")]  <- 1
		}
		for(i in 1:(n_disorders-1)){
		for(j in (i+1):n_disorders){
			temp_prsprs  <- crosstrait_cor.prs[i,j] * sqrt(prs_r2l[i]*prs_r2l[j])
			temp_prsliab <- cov_prsAliabB_theory(rg_AB=crosstrait_rg[i,j],  prs_r2l_A=prs_r2l[i] , h2l_A=snp_h2l[i] ,  h2l_B=snp_h2l[j])
			temp_liabprs <- cov_prsAliabB_theory(rg_AB=crosstrait_rg[i,j],  prs_r2l_A=prs_r2l[j] , h2l_A=snp_h2l[j] ,  h2l_B=snp_h2l[i])
			VARCOV[paste("prs_",i,sep="") ,paste("prs_",j,sep="") ] <- VARCOV[paste("prs_",j,sep="") ,paste("prs_",i,sep="") ] <- temp_prsprs		
			VARCOV[paste("prs_",i,sep="") ,paste("liab_",j,sep="")] <- VARCOV[paste("liab_",j,sep=""),paste("prs_",i,sep="") ] <- temp_prsliab		
			VARCOV[paste("liab_",i,sep=""),paste("prs_",j,sep="") ] <- VARCOV[paste("prs_",j,sep="") ,paste("liab_",i,sep="")] <- temp_liabprs		
			VARCOV[paste("liab_",i,sep=""),paste("liab_",j,sep="")] <- VARCOV[paste("liab_",j,sep=""),paste("liab_",i,sep="")] <- varcov_liab[i,j]	
		}}
		MEAN <- rep(0,each=dim(VARCOV)[1]) 

	#####################################
	## expected VARCOV & MEAN of prs and liability in truncated sub-populations
	#####################################

		truncate_mvnorm <- function( mu , sigma , i_truncate , threshold ){
			## mu: vector of means
			## sigma: dataframe of covariance
			## i_truncate: index variable to be truncated
			## threshold: threshold for truncating (var>threshold)
			## when wishing to select control, simply set threshold to be -1*threshold of cases

			## step 1: scale index variable to N(0,1)
			var_scaler       <- sqrt(1/sigma[i_truncate,i_truncate])
			mu_shift         <- var_scaler * mu
			sigma_scaled     <- var_scaler*var_scaler*sigma
			mu_scaled        <- rep(0,each=length(mu))
			threshold_scaled <- (threshold*var_scaler)-mu_shift[i_truncate]*sign(threshold) ## times sign(threshold) to allow for defining controls by negative thresholds

			## step 2: get updated {sigma_scaled & mu_scaled} based on truncating index variable
			K <- pnorm(threshold_scaled,lower.tail=FALSE)
			z <- dnorm(threshold_scaled) ## height_normal at threshold, following notation work with Naomi
			i <- z/K                     ## mean liability in cases, following notation work with Naomi
			k <- i*(i-threshold_scaled)  ## variance reduction in cases, following notation work with Naomi. Bulmer. The mathematical theory of quantitative genetics. 1985
			if(threshold<0){i <- -1*i} ## i.e. when selecting controls (note this does not impact k as above)
			help_matrix <- array( sigma_scaled[,i_truncate] , dim=dim(sigma_scaled))
			sigma_scaled_truncated <- sigma_scaled - k * help_matrix*t(help_matrix) ## Note: Tallis et al. 1987 Theor. Appl. Genet.: cov(X,Y|Z truncated)=cov(X,Y|Z not yet truncated)-k*cov(X,Z)*cov(Y,Z)  
			mu_scaled_truncated <- i * sigma_scaled[,i_truncate]/sigma_scaled[i_truncate,i_truncate] ## Note: cov/var is slop of regressing var on index_var

			## step 3: rescale to original scale
			sigma_truncated <- sigma_scaled_truncated/(var_scaler*var_scaler)
			mu_truncated <- (mu_scaled_truncated+mu_shift)/var_scaler

			return(list(sigma_truncated=sigma_truncated,mu_truncated=mu_truncated))
		} ## end truncate_mvnorm()

		mvnorm_list  <- list()
		mvnorm_list[[length(mvnorm_list)+1]] <- list(mu=MEAN,sigma=VARCOV) ; names(mvnorm_list)[[length(mvnorm_list)]] <- "population"
		for(row in 1:dim(liab.configuration)[1]){
			mu <- MEAN 
			sigma <- VARCOV
			for(i_disorder in 1:n_disorders){
				if(liab.configuration[row,analytical_col[i_disorder]]==1){
					temp <- truncate_mvnorm( mu=mu , sigma=sigma , i_truncate=i_disorder+length(analytical_col) , threshold=thresholds[i_disorder] ) ## prs preceed liability
					mu <- temp$mu_truncated
					sigma <- temp$sigma_truncated
				}
				if(liab.configuration[row,analytical_col[i_disorder]]==0){
					temp <- truncate_mvnorm( mu=mu , sigma=sigma , i_truncate=i_disorder+length(analytical_col) , threshold=-thresholds[i_disorder] ) ## prs preceed liability
					mu <- temp$mu_truncated
					sigma <- temp$sigma_truncated
				}
			}
			mvnorm_list[[length(mvnorm_list)+1]] <- list(mu=mu,sigma=sigma) ; names(mvnorm_list)[[length(mvnorm_list)]] <- liab.configuration$quick_name[row]
		}

	#####################################
	## Compute posterior probability
	#####################################

		f <- function(prs_liab_i){
			prs_liab_i <- as.numeric(prs_liab_i)
			for(row in 1:dim(liab.configuration)[1]){ w <- mvnorm_list[[liab.configuration$quick_name[row]]]; liab.configuration$dmvnorm[row] <- dmvn(X=prs_liab_i, mu=w$mu[1:length(analytical_col)], sigma=w$sigma[1:length(analytical_col),1:length(analytical_col)] )} 
			liab.configuration$test_priorprob*liab.configuration$dmvnorm/sum(liab.configuration$test_priorprob*liab.configuration$dmvnorm)
		}
		post_prob <- as.data.frame(t(apply(X=prs_liab,MAR=1,FUN=f))) ; colnames(post_prob) <- liab.configuration$quick_name
		for(category_name in category_names){
			quick_names <- liab.configuration$quick_name[ liab.configuration$diagnostic.category==category_name ]
			post_prob[ , category_name ] <- rowSums(post_prob[,quick_names,drop=FALSE])
		}
		post_prob <- post_prob[,c(category_names,liab.configuration$quick_name)]
		colnames(post_prob) <- paste("prob_",colnames(post_prob),sep="") 

		return(list(post_prob=as.data.frame(post_prob),liab.configuration=liab.configuration,mvnorm_list=mvnorm_list))

} ## end DDxPRS()

