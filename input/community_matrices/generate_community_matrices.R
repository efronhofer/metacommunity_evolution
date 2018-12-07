#####################################################################################
# generation of communty matrices
#
# Emanuel A. Fronhofer
# March 2018
#
####################################################################################
#
# Copyright (C) 2018  Emanuel A. Fronhofer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################
rm(list=ls())

####################################################################################
# GENERAL MODEL PARAMETERS
# no. species
n_species <- 80
# no. replicates
n_replicates <- 10
# structure of the matrix: only competitive inter-specific interactions ("competition") or allow also for facilitation ("all")
interaction_mode <- "competition"
# symmetry of the interaction matrix ("symmetric", "non-symmetric")
interaction_symmetry <- "symmetric"

####################################################################################
# SPECIES SPECIFIC PARAMETERS
# mean competition coefficient (alpha)
alpha_mean_interspecific <- 0.001
alpha_mean_intraspecific <- 0.002
# variance in competition coefficients
alpha_var_interspecific <- 0.1
alpha_var_intraspecific <- 0.1

# fraction of facilitative interactions
frac_facilitation <- 0

for (act_replicate in 1:n_replicates){
  # set up matrix
  # symmetric matrix
  if(interaction_symmetry == "symmetric"){
    alpha <- matrix(0,nrow=n_species, ncol=n_species)
    # make the interaction matrix alpha[effect on species Y, effect of species X]
    # generate the lower part and then add the transposed matrix
    if (interaction_mode == "all"){
      alpha[lower.tri(alpha)] <- rlnorm(length(alpha[lower.tri(alpha)]),log(alpha_mean_interspecific), sqrt(alpha_var_interspecific)) * sample(c(-1,1),length(alpha[lower.tri(alpha)]),replace=T,prob=c(frac_facilitation,1-frac_facilitation))
    }else{
      if(interaction_mode=="competition"){
        alpha[lower.tri(alpha)] <- rlnorm(length(alpha[lower.tri(alpha)]),log(alpha_mean_interspecific), sqrt(alpha_var_interspecific))
      }else{
        print("ERROR in interaction mode.")
      }
    }
    
    alpha <- alpha + t(alpha)
    
  }else{
    if(interaction_symmetry == "non-symmetric"){
      
      if (interaction_mode == "all"){
        alpha <- matrix(rlnorm(n_species^2,log(alpha_mean_interspecific), sqrt(alpha_var_interspecific)),nrow=n_species, ncol=n_species) * matrix(sample(c(-1,1),n_species^2,replace=T,prob=c(frac_facilitation,1-frac_facilitation)) ,nrow=n_species, ncol=n_species)
      }else{
        if(interaction_mode=="competition"){
          alpha <- matrix(rlnorm(n_species^2,log(alpha_mean_interspecific), sqrt(alpha_var_interspecific)),nrow=n_species, ncol=n_species)
        }else{
          print("ERROR in interaction mode.")
        }
      }
      
    }else{print("ERROR in matrix symmetry definition.")}
  }

  # generate the diagonal values separately
  # these values should be positive only so they are drawn from a lognormal distribution
  diag(alpha) <- rlnorm(length(diag(alpha)),log(alpha_mean_intraspecific), sqrt(alpha_var_intraspecific))

  write.table(alpha,file=paste("community_matrix_",act_replicate,"_",interaction_mode,"_",interaction_symmetry,".txt", sep=""), quote = F, row.names = F, col.names = F)
}
