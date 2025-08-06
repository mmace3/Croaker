
get_data <- function(obs, fit){
  if(max(obs) < max(fit))
  {
    fit_hist <- hist(fit, plot = FALSE)
    mids <- fit_hist$mids
    brks <- fit_hist$breaks
    fit_counts <- fit_hist$counts
    obs_hist <- hist(obs, breaks = brks, plot = FALSE)
    obs_counts <- obs_hist$counts
    length_brks <- length(brks)
    low_brks <- brks[-length_brks]
    high_brks <- brks[-1]
    
    res <- data.frame(mids = mids,
                      obs = obs_counts,
                      fit = fit_counts,
                      sqrt_obs = sqrt(obs_counts),
                      sqrt_fit = sqrt(fit_counts),
                      low_brks = low_brks,
                      high_brks = high_brks
    )
  } else {
    
    obs_hist <- hist(obs, plot = FALSE)
    mids <- obs_hist$mids
    brks <- obs_hist$breaks
    obs_counts <- obs_hist$counts
    fit_hist <- hist(fit, breaks = brks, plot = FALSE)
    fit_counts <- fit_hist$counts
    length_brks <- length(brks)
    low_brks <- brks[-1]
    high_brks <- brks[-length_brks]
    
    res <- data.frame(mids = mids,
                      obs = obs_counts,
                      fit = fit_counts,
                      sqrt_obs = sqrt(obs_counts),
                      sqrt_fit = sqrt(fit_counts),
                      low_brks = brks[-length_brks],
                      high_brks = brks[-1]
    )
  }
  
  return(res)
}