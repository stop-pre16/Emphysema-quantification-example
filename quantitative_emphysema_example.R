####   R script with example code for estimating various quantitative emphysema measures using a simulated 
####   3D point pattern

library(dplyr)
library(car)
library(pbapply)
library(sncp)
library(spatstat)
library(igraph)

###  Read in simulated point pattern provided
sim_pattern = readRDS(file = "local/path/to/this/file/sim_pattern_3d.RDS")

###   Visualize 3D pattern
car::scatter3d(x = sim_pattern[, 1], 
               y = sim_pattern[, 2],
               z =sim_pattern[, 3],
               surface = F, 
               point.col = "black")

###  This pixelated pattern was simulated in a cube with dimensions [0, 99] x [0, 99] x [0, 99]
ubound = 99

xwin = c(0, ubound)
ywin = c(0, ubound)
zwin = c(0, ubound)

###  First we'll calculate the power law exponent D from Mishima et al (1999) [1]

# Create 3D point pattern object for spatstat functions
dat_pp3 = pp3(x = sim_pattern$x, y = sim_pattern$y, z = sim_pattern$z, box3(xwin, ywin, zwin))

# Identify contiguous LAA cluster using a connected components function from spatstat. 
# We select R = 1.8 to allow for diagonal connections
res_cc = connected.pp3(X = dat_pp3, R = 1.8)

# Calculate the sizes in terms of number of voxels for all 3D clusters
cc_sum = res_cc$data %>% data.frame %>% group_by(marks) %>% summarize(npts = n()) %>% arrange(desc(npts)) %>% data.frame

# Estimate the power law exponent D using a function from igraph package and the MLE implementation
fmle = fit_power_law(x = cc_sum$npts, implementation = "R.mle")
estimated_D = fmle@coef
print(estimated_D)

################################################################################################
################################################################################################
################################################################################################
################################################################################################

###   Next we calculate Normalized Join-Count (NJC) from Virdee et al (2021) [2]

# Calculating total joins given our observed space is a 3D cube is straightforward

# Joins within each 2D slice:
joins_all_within_slice = (99*100) * 2

# Joins between 2 adjacent 2D slices: 
joins_all_between_slice = 100 * 100

# Total joins = joins within a slice times number of slices + joins between slices * number of slices minus 1
joins_all_total = joins_all_within_slice *100 + joins_all_between_slice *99

###  To calculate LAA-LAA joins, we go sequentially through the slices and compute within- and between-slice joins 
###  Note this same strategy can be also be used to compute total joins if the observed space is more irregular
###  like a 3D segmented lung
joins_laa_each_slice = bind_rows(pblapply(0:ubound, function(zz){
  
  dsub = sim_pattern %>% filter(z == zz) # Filter to LAA voxels within the given slice
  dsub2 = sim_pattern %>% filter(z == (zz+1)) # Filter to LAA voxels in the slice above the current one
  cnt1_e = 0 # Within-slice join counter
  cnt2_f = 0 # Between-slice join counter
  
  # Compute pairwise distance matrix for all LAA voxels in current slice, then sum up unique adjacencies
  # based on shared edges
  if(nrow(dsub) > 1){
    dmat1 = fastPdist2(A = as.matrix(dsub[, 1:2], ncol = 2), B = as.matrix(dsub[, 1:2], ncol = 2))
    dvec1 = dmat1[upper.tri(dmat1, diag = F)]
    cnt1_e = sum(dvec1 < 1.01)
  }
  # Compute pairwise distances between LAA voxels in current slice and LAA voxels in the next slice, then sum up
  # unique adjacencies based on shared faces
  if(nrow(dsub2) > 0){
    dmat2 = fastPdist2(A = as.matrix(dsub[, 1:3], ncol = 3), B = as.matrix(dsub2[, 1:3], ncol = 3))
    cnt2_f = sum(dmat2 < 1.01)
  }
  ret = data.frame(z = zz, 
                   cnt_f = cnt1_e + cnt2_f)
  return(ret)
}))

# Sum all LAA-LAA joins
joins_laa_total = sum(joins_laa_each_slice$cnt_f)

# Compute NJC as LAA-LAA joins divided by total joins
NJC = joins_laa_total / joins_all_total

################################################################################################
################################################################################################
################################################################################################
################################################################################################

###   Finally, we apply the spatial birth-death MCMC implemented in the sncp R package to compute the spatial point
###   process measures described in Vestal et al (2019) [3]

n_it = 10000

# Loop through each axial slice and apply a 2D model in each (will likely take 10+ minutes)
res_sncp_all_slices = bind_rows(pblapply(0:ubound, function(zz){
  dsub = sim_pattern %>% filter(z == zz)
  if(nrow(dsub) > 0){
    msub = cbind(as.numeric(dsub$x), as.numeric(dsub$y))
    
    # Apply the BD-MCMC using similar priors/parameters as described in [3]
    res_bdmcmc = sncp_bdmcmc_cont(obs_points = msub,
                                  obs_window = as.matrix(rbind(xwin, ywin)), 
                                  LM = ubound^2,
                                  mean_mu_alpha = 4.5, 
                                  sd_log_alpha = 1, 
                                  sd_prop_alpha = .5, 
                                  beta = .001, 
                                  n_it = n_it, 
                                  window_hw = 5, 
                                  df_iw_prior = 5, 
                                  df_iw_prop = 10, 
                                  sigma_prior = diag(10, nrow = 2), 
                                  var_mu_alpha = 3,
                                  pen_dist = 10, 
                                  pen_val = 1e-9, 
                                  n_cent_init = 3, 
                                  prior_n_cent = 10,
                                  max_bd_events = 15, 
                                  max_bd_vt = 50)
    
    # Discard the 1st 10% of the MCMC samples as burn-in
    bbound = round(n_it * .10) + 1
    
    # Estimate number of centers, average cluster size, and noise/diffuse intensity as the posterior medians
    est_nc = median(res_bdmcmc$n_centers_sample[bbound:n_it])
    est_ma = median(res_bdmcmc$mu_alpha_sample[bbound:n_it])
    est_beta =  median(res_bdmcmc$beta_sample[bbound:n_it])
    res = data.frame(z = zz, 
                     est_nc = est_nc, 
                     est_ma = est_ma, 
                     est_beta = est_beta)
  }
  else{
    res = data.frame(z = zz, 
                     est_nc = NA, 
                     est_ma = NA, 
                     est_beta = NA)
  }
  return(res)
}))

head(res_sncp_all_slices)

# Generate overall summary as the mean across all 2D slices
sncp_summary = data.frame(mean_NC = mean(res_sncp_all_slices$est_nc, na.rm = T),
                          mean_ACS = mean(res_sncp_all_slices$est_ma, na.rm = T),
                          mean_diffuse = mean(res_sncp_all_slices$est_beta, na.rm = T))
sncp_summary

################################################################################################
################################################################################################
################################################################################################
################################################################################################

###   References
# [1] Mishima M, Hirai T, Itoh H, Nakano Y, Sakai H, Muro S, Nishimura K, Oku Y, Chin K, Ohi M, Nakamura T. 
# Complexity of terminal airspace geometry assessed by lung computed tomography in normal subjects and patients with chronic obstructive pulmonary disease.
# Proceedings of the National Academy of Sciences. 1999 Aug 3;96(16):8829-34.

# [2] Virdee S, Tan WC, Hogg JC, Bourbeau J, Hague CJ, Leipsic JA, Kirby M. 
# Spatial dependence of CT emphysema in chronic obstructive pulmonary disease quantified by using join-count statistics.
# Radiology. 2021 Dec;301(3):702-9.

# [3] Vestal BE, Carlson NE, Estépar RS, Fingerlin T, Ghosh D, Kechris K, Lynch D. 
# Using a spatial point process framework to characterize lung computed tomography scans. 
# Spatial statistics. 2019 Mar 1;29:243-67.




