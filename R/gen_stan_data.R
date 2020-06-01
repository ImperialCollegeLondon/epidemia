

#' Generates data to pass to rstan::sampling
#' 
#' @param formula An R object of class `formula`. The left hand side must take the form `Rt(group,date)', with 'group' representing a factor vector indicating group membership (i.e. country, state, age cohort), and 'code' being a vector of Date objects.
#' @param data A dataframe with columns corresponding to the terms appearing in 'formula'. See [lm].
#' @param obs A named list giving available observations
#' * deaths: A three column dataframe representing death data. The first column represents group membership and must be coercible to class 'factor'. The second column indicates the observation date and must be coercible to class `Date'.
#' * incidence: Same as 'deaths', although giving incidence data.
#' * dtd: A vector representing 'days to observation'. Like si, this is a probability vector. The nth element giving the probability that the observation event (incidence/death) occurs n days after an individual is infected.
#' * dti: same as 'dtd', but representing time until an incidence is recorded after onset of infection.
#' @param pops  A two column dataframe giving the total population of each group. First column represents the group, while the second the corresponding population.
#' @param ifr A two column dataframe giving the infection fatality rate in each group. First column represents the group, while the second the corresponding IFR.
#' @param si A vector representing the serial interval of the disease (a probability vector).
#' @param seed_days Number of days for which to seed infections.
#' @param ... Arguments allowed in rstanarm::stan_glmer(). For example one can control the prior distribution of the covariates.
#' @examples
#' @return A list with required data to pass to rstan::sampling.
genStanData <- 
  function(formula, 
           data = NULL,
           obs,
           pops,
           ifr,
           si,
           seed_days = 6,
           ...) {

  dots     <- list(...)

  formula <- checkFormula(formula)
  data    <- checkData(formula, data)

  groups <- levels(data$group)

  obs   <- checkObs(data, obs)
  pops  <- checkPops(pops, groups)
  ifr   <- checkIFR(ifr, groups)
  si    <- checkSV(si)

  if (seed_days < 1)
    stop("'seed_days' must be greater than zero", call. = FALSE)

  # standata <- genModelStanData(levels(covariates$group), 
  #                              data, 
  #                              obs_type, 
  #                              seed_days, 
  #                              forecast
  # )
  
  # # use only covariates within the simulation range
  # obs           <- standata$obs
  # f             <- function(x) cbind(x, idx = seq.int(nrow(x)))
  # obs           <- data.table::rbindlist(Map(f, split(obs, obs$group)))
  # covariates    <- dplyr::left_join(obs, covariates, by = c("group", "date"))
  # covariates    <- covariates[complete.cases(covariates),]
  # standata$memb <- as.numeric(covariates$group)
  # standata$idx  <- covariates$idx

  # # get all stan data relating to data$covariates
  # standata <- c(standata, 
  #               genCovariatesStanData(formula, 
  #                                     covariates, 
  #                                     ...
  #               )
  # )
                                
  return(nlist(data, obs))
} 




















# #' Generates data to pass to rstan::sampling
# #' 
# #' @param formula An R object of class `formula`. The left hand side must take the form `cbind(y1,y2)', with 'y1' representing a factor vector indicating group membership (i.e. country, state, age cohort), and 'y2' being a vector of Date objects.
# #' @param data A list describing all data required for fitting 
# #' * obs: A three column dataframe representing the response data (i.e. incidence or death counts) . The first column represents group membership and must be coercible to class 'factor'. The second column indicates the observation date and must be coercible to class `Date'. The final column is the observed data itself.
# #' * pops: A two column dataframe giving the total population of each group. Again, first column represents the group, while the second the corresponding population.
# #' * ifr: A two column dataframe giving the infection fatality rate in each group. First column represents the group, while the second the corresponding IFR.
# #' * si: A vector representing the serial interval of the disease (a probability vector).
# #' * dto: A vector representing 'days to observation'. Like si, this is a probability vector. The nth element giving the probability that the observation event (incidence/death) occurs n days after an individual is infected.
# #' * covariates: An optional (but recommended) list element containing a dataframe. This dataframe must include a column corresponding to each term in 'formula'. If not included, then the terms will be looked for in the environment corresponding to 'formula'.
# #' @param obs_type Currently either "Incidence" or "Deaths" are supported.
# #' @param seed_days Number of days for which to seed infections.
# #' @param forecast Number of days to continue the simulation after the final date in data$obs 
# #' @param ... Arguments allowed in rstanarm::stan_glmer(). For example one can control the prior distribution of the covariates.
# #' @examples
# #' @return A list with required data to pass to rstan::sampling.
# genStanData <- 
#   function(formula, 
#            data, 
#            obs_type = c("Incidence", "Deaths"),  
#            seed_days = 6, 
#            forecast = 7, 
#            ...) {

#   obs_type <- match.arg(obs_type)
#   dots     <- list(...)

#   formula <- as.formula(formula)

#   lhs <- formula[[2]]
#   # y   <- getResponses(lhs, data)
#   # Assuming for now that data$covariates is not null
#   # Todo: allow data$covariates to be null
#   covariates <- processCovariates(lhs, data$covariates)

#   # check and manipulate data argument
#   data <- checkData(data, levels(covariates$group))

#   if (!length(obs_type))
#     stop("'obs_type' must be one of ", paste(obs_types, collapse = ", "))
#   if (seed_days < 1)
#     stop("'seed_days' must be greater than zero")
#   if (forecast < 1)
#     stop("'forecast' must be greater than zero")

#   standata <- genModelStanData(levels(covariates$group), 
#                                data, 
#                                obs_type, 
#                                seed_days, 
#                                forecast
#   )
  
#   # use only covariates within the simulation range
#   obs           <- standata$obs
#   f             <- function(x) cbind(x, idx = seq.int(nrow(x)))
#   obs           <- data.table::rbindlist(Map(f, split(obs, obs$group)))
#   covariates    <- dplyr::left_join(obs, covariates, by = c("group", "date"))
#   covariates    <- covariates[complete.cases(covariates),]
#   standata$memb <- as.numeric(covariates$group)
#   standata$idx  <- covariates$idx

#   # get all stan data relating to data$covariates
#   standata <- c(standata, 
#                 genCovariatesStanData(formula, 
#                                       covariates, 
#                                       ...
#                 )
#   )
                                
#   return(standata)
# } 


# processCovariates <- function(lhs, covariates) {
#   if (lhs[[1]] != "cbind")
#     stop("Left hand side of 'formula' must be of the form 'cbind(y1,y2)'")

#   covariates$group <- covariates[[lhs[[2]]]]
#   covariates$dates  <- covariates[[lhs[[3]]]]

#   if(is.null(covariates$group)) 
#     stop(paste0("cannot find ", lhs[[2]], "in ", env, call.=FALSE))
#   if(is.null(covariates$dates)) 
#     stop(paste0("cannot find ", lhs[[3]], "in ", env, call.=FALSE))
  
#   covariates$group  <- as.factor(covariates$group)
#   covariates$dates   <- as.Date(covariates$dates)

#   return(covariates)
# }



# Get the two vectors referred to in the LHS of 'formula' argument of genStanData
#
# @param lhs The left hand side of 'formula' object.
# # @param data Named list. See [genStanData]
# getResponses <- function(lhs, data) {
#   if (lhs[[1]] != "cbind")
#     stop("Left hand side of 'formula' must be of the form 'cbind(y1,y2)'")

#   env <- if (is.null(data$covariates))  environment(formula) else  data$covariates

#   groups <- env[[lhs[[2]]]]
#   dates  <- env[[lhs[[3]]]]

#   if(is.null(groups)) 
#     stop(paste0("cannot find ", lhs[[2]], "in ", env, call.=FALSE))
#   if(is.null(dates)) 
#     stop(paste0("cannot find ", lhs[[3]], "in ", env, call.=FALSE))
  
#   groups  <- as.factor(groups)
#   dates   <- as.Date(dates)

#   return(nlist(groups, dates))
# }



