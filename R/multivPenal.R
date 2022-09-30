#' Competing Joint Frailty Model: A single type of recurrent event and two
#' types of terminal events.
#'
#' @description Fit a joint frailty model for a single recurrent event
#' and two terminal events using a parametric weibull model.
#'
#' @aliases multivPenal transfo.table multivPenal for multivariate frailty
#' model
#'
#' @usage
#' multivPenal(
#' formula,
#' formula.terminalEvent,
#' formula.terminalEvent2,
#' data,
#' initialize = TRUE,
#' maxit = 350,
#' gapTimes=FALSE,
#' jointGeneral=FALSE
#' hazard = "Weibull",
#' GHpoints = 32,
#' tolerance = rep(10^-3, 3),
#' init.hazard = c(1,0.5,1,1,1,1),
#' init.Sigma = 0.5,
#' init.Alpha1 = 0.1,
#' init.Alpha2 = -0.1,
#' init.B = c(0,0,0))
#'
#' @param formula a formula object, with the response for the first recurrent
#' event on the left of a \eqn{\sim} operator, and the terms on the right.
#' The response must be in the format Surv(t0, t1, recurrentevent),
#' where t0 is the start time for an at-risk period for the recurrent event,
#' t1 is the end time for an at-risk period for the recurrent event, and
#' recurrent event is a numeric indicator for whether an event was observed (1)
#' or was censored (2).
#' This formulation is used regardless of whether `gapTimes` is true or false.
#' If gapTimes is true, the duration of the risk period will be calculated as
#' t1 - t0.
#'
#' @param formula.terminalEvent, a formula object,
#' empty on the left of a \eqn{\sim} operator,
#' and the terms on the right. Leave the formula at
#' the default value (NULL) for a model with no variables.
#'
#' @param formula.terminalEvent2, a formula object,
#' empty on the left of a \eqn{\sim} operator,
#' and the terms on the right.Leave the formula at
#' the default value (NULL) for a model with no variables.
#'
#' @param data a 'data.frame' with the variables used in 'formula',
#' 'formula.terminalEvent', and 'formula.terminalEvent2'.
#'
#' @param initialize Logical value to internally initialize regression coefficients and
#' baseline hazard functions parameters using simpler models from frailtypack.
#' When initialization is requested, the
#' program first fits two joint frailty models for the recurrent
#' events and each terminal event.
#' When FALSE, parameters are initialized via the arguments
#' init.hazard, init.Sigma, init.Alpha1, init.Alpha2, init.B.
#'
#' @param gapTimes Should the hazard be formulated at time from the beginning of the
#' risk period? If FALSE, time starts at t0 = 0.
#'
#' @param jointGeneral Logical. Should the model have separate (correlated)
#' frailty random effects that connect each terminal event to the recurrent
#' event? Default is FALSE.
#'
#' @param maxit maximum number of iterations for the Marquardt algorithm.
#' Default is 350.
#'
#' @param hazard Type of hazard functions: For now the only option is "Weibull"
#' for parametric Weibull function.
#'
#' @param GHpoints Integer. Number of nodes for Gauss-Hermite integration
#' to marginalize random effects/frailties. Default is 32.
#'
#' @param tolerance Numeric, length 3. Optimizer's tolerance for (1) successive change
#' in parameter values, (2) log likelihood, and (3) score, respectively.
#'
#' @param init.hazard Numeric. Initialization values for hazard parameters.
#' If a weibull model is used, the order is:
#' shapeR, scaleR, shapeTerminal1, scaleTerminal1, shapeTerminal2, scaleTerminal2.
#'
#' @param init.Sigma Numeric,. Initialization value for the standard deviation of the
#' normally-distributed random effects.
#'
#' @param init.Alpha1 Numeric. Initialization value for the parameter alpha that
#' links the hazard function of the recurrent event to the first terminal event.
#'
#' @param init.Alpha2 Numeric. Initialization value for the parameter alpha that
#' links the hazard function of the recurrent event to the second terminal event.
#'
#' @param init.B Numeric vector, same length and order as the covariate vectors
#' for each the recurrent event, terminal1, and terminal2 (in that order).
#'
#' @details{
#' \if{html}{Right-censored data are allowed.
#' Left-truncated data and stratified analysis are not possible.}
#' }
#'
#' @return Parameters estimates of a competing joint frailty model, more
#' generally a 'competingPenal' object. Methods defined for 'competingPenal' objects
#' are provided for print, plot and summary. The following components are
#' included in a 'competingPenal' object for multivariate Joint frailty models.
#'
#' \item{summary.table}{A table describing the estimate, standard error,
#' confidence interval, and pvalues for each of the parameters in the model.}
#'
#' \item{b}{sequence of the corresponding estimation of the splines
#' coefficients, the random effects variances, the coefficients of the
#' frailties and the regression coefficients.}
#'
#' \item{call}{The code used for fitting the model.}
#'
#' \item{n}{the number of observations used in the fit.}
#'
#' \item{groups}{the number of subjects used in the fit.}
#'
#' \item{n.events}{the number of recurrent events of type 1 observed in the fit.}
#'
#' \item{n.terminal1}{the number of terminal1 events observed in the fit.}
#'
#' \item{n.terminal2}{the number of terminal2 events observed in the fit.}
#'
#' \item{loglik}{the marginal log-likelihood in the parametric case.}
#'
#' \item{AIC}{the Akaike information Criterion for the parametric
#' case.\deqn{AIC=\frac{1}{n}(np - l(.))}}
#'
#' \item{npar}{total number of parameters.}
#'
#' \item{nvar}{A vector with the number of covariates of each type of hazard function as
#' components.}
#'
#' \item{varH}{the variance matrix of all parameters before
#' positivity constraint transformation (theta and
#' the hazard coefficients).
#' Then, the delta method is needed to obtain the estimated variance parameters.}
#'
#' \item{varHIH}{the robust estimation of the
#' variance matrix of all parameters (theta and
#' the hazard coefficients).}
#'
#' \item{formula}{the formula part of the code used
#' for the model for the recurrent event.}
#'
#' \item{formula.terminalEvent}{the formula part of the code used for the model
#' for the terminal1 event.}
#'
#' \item{formula.terminalEvent2}{the formula part of the code used for the model
#' for the terminal2 event.}
#'
#' \item{x1}{vector of times for hazard functions of
#' the recurrent events are estimated.
#' By default seq(0,max(time),length=99),
#' where time is the vector of survival times.}
#'
#' \item{lam1}{matrix of hazard estimates and confidence bands for recurrent
#' events of type 1.}
#'
#' \item{xSu1}{vector of times for the survival function of
#' the recurrent event of type 1.}
#'
#' \item{surv1}{matrix of baseline survival
#' estimates and confidence bands for recurrent events of type 1.}
#'
#' \item{x3}{Vector of times for terminal1 (see x1 value).}
#'
#' \item{lam3}{Matrix of hazard, CI for terminal1 evaluated at times x3.}
#'
#' \item{xSu3}{Vector of times for the survival function of terminal1.}
#'
#' \item{surv3}{Vector of the survival function of terminal1 evaluated at xSu3.}
#'
#' \item{x4}{Vector of times for terminal2 (see x1 value).}
#'
#' \item{lam4}{Matrix of hazard, CI for terminal1 evaluated at times x2.}
#'
#' \item{xSu4}{Vector of times for the survival function of terminal2.}
#'
#' \item{surv4}{Vector of the survival function of terminal1 evaluated at xSu2.}
#'
#' \item{median1}{The value of the median survival and its confidence bands for the recurrent event of type 1.}
#'
#' \item{median3}{The value of the median survival and its confidence bands for terminal1.}
#'
#' \item{median4}{The value of the median survival and its confidence bands for the terminal2.}
#'
#' \item{n.iter}{number of iterations needed to converge.}
#'
#' \item{noVar}{indicator vector for recurrent, death and recurrent 2
#' explanatory variables.}
#'
#' \item{istop}{Vector of the convergence criteria.}
#'
#' \item{martingale.res}{martingale residuals for each cluster (recurrent).}
#'
#' \item{martingale2.res}{martingale residuals for each cluster
#' (recurrent of type 2).}
#'
#' \item{martingale4.res}{martingale residuals for
#' each cluster (death).}
#'
#' \item{frailty.pred}{empirical Bayes prediction of the
#' frailty term(s).}
#'
#' \item{frailty.var}{variance of the empirical Bayes
#' prediction of the frailty term(s).}
#'
#' \item{frailty.corr}{Correlation between the empirical Bayes prediction of
#' the two frailty.}
#'
#' \item{linear.pred}{Predicted log relative-hazard of the recurrent event.}
#'
#' \item{linear.pred3}{Predicted log relative-hazard of the terminal1 event.}
#'
#' \item{linear.pred4}{Predicted log relative-hazard of the terminal2 event.}
#'
#' @seealso
#' \code{\link{terminal}},
#' \code{\link{terminal2}},
#'
#' @export
#'
#'
"multivPenal" <-
  function(formula,
           formula.terminalEvent = NULL,
           formula.terminalEvent2 = NULL,
           data,
           initialize = TRUE,
           gapTimes = FALSE,
           jointGeneral = FALSE,
  		   maxit = 350,
           hazard = "Weibull",
           GHpoints = 32,
           tolerance = rep(10^-3, 3),
           init.hazard = NULL,
           init.Sigma = 0.5,
           init.Alpha1 = 0.1,
           init.Alpha2 = -0.1,
           init.B = NULL)
{
#############################################################
#############################################################
# Outline of competingPenal.R
# (1) Verify/extract event indicators, times, groups
# (2) Configure Hazards
# (3) Configure Model Matrices
# (4) Configure Parameters
# (5) Define Gauss-Hermite Nodes and Weights
# (6) Compute Gap times, if gapTimes == TRUE
# (7) Fill Starting Parameter Vector with User-Defined Values
#		OR Initialize Models
# (8) Check Dimensions of all variables being sent to Fortran
#       for debugging
# (9) Send to Fortran for optimization
# (10) Format Model Summary Table
# (11) Format Initialization Model Summary Tables
#		(if initialization == TRUE)
#############################################################
#############################################################

n.knots <- 0 # remove when splines models are available
kappa <- rep(0, 4)
crossVal <- 0

if (class(formula) != "formula")
  stop("The argument formula must be a formula")

if (nrow(data) == 0)
  stop("No (non-missing) observations")

#############################################################
#############################################################
# (1) Verify/extract event indicators, times, groups

### (1a) Verify Recurrent Event 1 / Times
# NN =  names of time gap and recurrent event 1
Y1 <- get_all_vars(update(formula, "~1"), data)

if (ncol(Y1) != 3) {
  stop(
  	"Survival object outcome must be specified using \n
  	time when at risk period starts and stops, such as: \n
  	Surv(time.start, time.stop, event.indicator).
  	This is true for both calendar time and gap time
  	formulations."
  )
}

NN <- colnames(Y1)

EVENT1 <- NN[3]
TSTART <- NN[1]
TSTOP <- NN[2]

event1 <- Y1[, 3]
tt10 <- Y1[, 1]
tt11 <- Y1[, 2]

### (1b) Verify recurrent event 2 (Not yet implemented)
event2.ind <- 0
event2 <- 0
tt0meta0 <- 0
tt1meta0 <- 0

# if (!missing(formula.Event2)) {
#   Y2 <- get_all_vars(update(formula.Event2, "~1"), data.Event2)
#
#   if (ncol(Y2) != 3) {
#   	stop(
#   		"Survival object outcome must be specified using \n
#   		time when at risk period starts and stops, such as: \n
#   		Surv(time.start, time.stop, event.indicator).
#   		This is true for both calendar time and gap time
#   		formulations."
#   	)
#   }
#   NN <- colnames(Y2)
#   EVENT2 <- NN[3]
#   TSTART2 <- NN[1]
#   TSTOP2 <- NN[2]
#
#   event2 <- Y2[, 3]
#   tt0meta0 <- Y2[, 1]
#   tt1meta0 <- Y2[, 2]
#
#   if (!all(event2 %in% c(1, 0))) {
#   	stop(
#   		"event2 must contain a variable coded 0-1 and a non-factor variable"
#   	)
#   }
#   event2.ind <- 1
# } else{
#   event2.ind <- 0
#   event2 <- 0
#   tt0meta0 <- 0
#   tt1meta0 <- 0
# }

### (1c) Verify Terminal Event 1
TT <-
  survival::untangle.specials(terms(formula, c("terminal")), "terminal", 1:10)$vars
start <- grep("\\(", unlist(strsplit(TT, "")))
stop <- grep("\\)", unlist(strsplit(TT, "")))
TERMINAL1 <- substr(TT, start = start + 1, stop = stop - 1)
if (length(TERMINAL1) == 0) {
  stop("A term for a terminal event must be included in the formula.")
}
if (!all(data[[TERMINAL1]] %in% c(1, 0))) {
  stop("terminal must contain a variable coded 0-1 and a non-factor variable")
}

### (1d) Verify Terminal Event 2
TT <-survival::untangle.specials(terms(formula, c("terminal2")), "terminal2", 1:10)$vars
start <- grep("\\(", unlist(strsplit(TT, "")))
stop <- grep("\\)", unlist(strsplit(TT, "")))
TERMINAL2 <- substr(TT, start = start + 1, stop = stop - 1)
if (length(TERMINAL2) == 0) {
  terminal2.ind <- 0
} else{
  if (!all(data[[TERMINAL2]] %in% c(1, 0))) {
  	stop(
  		"terminal must contain a variable coded 0-1 and a non-factor variable"
  	)
  }
  terminal2.ind <- 1
}

### (1e) Verify Cluster Variable
# name of cluster variable
TT <-
  survival::untangle.specials(terms(formula, c("cluster")), "cluster", 1:10)$vars
start <- grep("\\(", unlist(strsplit(TT, "")))
stop <- grep("\\)", unlist(strsplit(TT, "")))
CLUSTER <- substr(TT, start = start + 1, stop = stop - 1)

if (length(data[[CLUSTER]])) {
  uni.cluster <- unique(data[[CLUSTER]])
} else{
  stop("grouping variable is needed")
}
if (length(uni.cluster) == 1) {
  stop("grouping variable must have more than 1 level")
}
if (event2.ind == 1){
	if(CLUSTER %in% colnames(data.Event2)) {
		stop("grouping variable must be present in data.Event2")
	}
	if (!all(uni.cluster %in% data.Event2[[CLUSTER]])) {
		stop("all groups must be represented in data.Event2")
	}
}
##########################################################################
##########################################################################
# (2) Configure Hazards

if (hazard != "Weibull") {
  stop("Only 'Weibull' has been implemented in competing joint model.")
}

haztemp <- hazard
#if (!hazard %in% c("Weibull", "Splines")) {
#  stop("Only 'Weibull' or 'Splines' hazard can be specified in hazard argument.")
#}

# typeof is the numerical indicator for hazard
typeof <- switch(hazard,
  	     "Splines" = 0,
  	     "Weibull" = 2)

# Specify whether spline points are regularly spaced or at percentiles
if(typeof == 0) {
  	equidistant <- 1
}else if(typeof == 2){
  	### Weibull
  	equidistant <- 1
}

#### Configure Splines Hazard (knots, pentalty), Not yet implemented
# if (typeof == 0) {
#   crossVal <- 1 # always do cross validation for splines
#   if (missing(kappa))
#   	stop("smoothing parameter (kappa1) is required")
#   if (missing(n.knots))
#   	stop("number of knots are required")
#   if (class(kappa) != "numeric")
#   	stop("The argument kappa must be a numeric")
#   if (class(n.knots) != "integer"){
#   	warning("Converting n.knots to integers")
#   	n.knots <- as.integer(n.knots)
#   }
#   if (length(n.knots) != 1) {
#   	stop("length of knots must be 1.")
#   }
#   if (length(kappa) != 2 + event2.ind + terminal2.ind) {
#   	stop(
#   		"length of kappa must be equal to the number of formulas
#   		for different event types (3 or 4) in the order
#   		recurrent1, recurrent2, terminal1, terminal2."
#   	)
#   }
#   if (event2.ind == 1 & terminal2.ind == 0) {
#   	kappa = c(kappa0, 0)
#   }
#   if (event2.ind == 0 & terminal2.ind == 1) {
#   	kappa = c(kappa[1:2], 0, kappa[3])
#   }
#   n.knots[n.knots < 4 & n.knots != 0] <- 4
#   n.knots[n.knots > 20] <- 20
# }else if(typeof == 2){
#   # if the hazard is weibull
#   if (!(missing(n.knots)) || !(missing(kappa))) {
#   	warning("When parametric hazard is not 'Splines'
#   		'kappa' and 'n.knots' arguments are ignored.")
#   }
#   n.knots <- 0
#   kappa <- rep(0, 4)
#   crossVal <- 0
# }

# End Hazard Configuration
#########################################################################
#########################################################################
# (3) Configure Model Matrices

# noVarEvent indicates whether there are no explanatory variables for the
# recurrent1, terminal1, recurrent2, and terminal2 events (in that order)
noVarEvent = c(0,0,1,0)

# delete specials from the formula to get just covariates
# on right side.
specials = c("strata", "cluster", "terminal", "event2", "terminal2")
Terms = terms(formula, specials = specials)

# Check if all terms on right side of Recurrent event formula are specials
if(length(unlist(attr(Terms, "specials"))) == (length(unlist(attr(Terms, "variables")))-2)){
	noVarEvent[1] <- 1
}

# Check for terminal1 formula
if(is.null(formula.terminalEvent)){
	noVarEvent[2] <- 1
}

# Check for terminal2 formula
if(is.null(formula.terminalEvent2)){
	noVarEvent[4] <- 1
}

# Recurrent Event 1 Model Matrix
if(noVarEvent[1] == 0){
modelmatrix1 =
  model.matrix(update(
  	drop.terms(
  		termobj = terms(formula),
  		unlist(attr(Terms, "specials")) - 1,
  		keep.response = TRUE
  	),
  	~ . - 1
  ),
  data)
}else{
	  modelmatrix1 = matrix(0)
}

# need to get densely-ranked ids
group1 <- as.numeric(factor(data[[CLUSTER]]))

# Compute Event Counts
# nevents1 <- tapply(event1, group1, sum)

# Recurrent Event 2 Model Matrix
if (event2.ind == 1) {
  group2 <- as.numeric(factor(data.Event2[[CLUSTER]]))
  # Compute Event Counts
  #nevents2 <- tapply(event2, group2, sum)
  Terms2 = terms(formula.Event2, specials = specials)
  if (!is.null(unlist(attr(Terms2, "specials")))) {
  	modelmatrix3 =
  		model.matrix(update(
  			drop.terms(
  				termobj = terms(formula.Event2),
  				unlist(attr(
  					Terms2, "specials"
  				)) - 1,
  				keep.response = TRUE
  			),
  			~ . - 1
  		),
  		data.Event2)
  } else{
  	modelmatrix3 =
  		model.matrix(update(formula.Event2, ~ . - 1),
  			 data.Event2)
  }
} else{
  modelmatrix3 = matrix(0)
  group2 = 0
}

# Terminal Event 1 Model Matrix
data.terminal <- do.call(what = "rbind",
  	 lapply(split(x = data, f = data[[CLUSTER]]),
  	        function(df) {
  	        	subset(df, df[[TSTOP]] == max(df[[TSTOP]]))
  	        }))
groupdc <- as.numeric(factor(data.terminal[[CLUSTER]]))
tt1dc <- data.terminal[[TSTOP]]
terminal1 <- data.terminal[[TERMINAL1]]

if(noVarEvent[2]==0){
modelmatrix2 = model.matrix(update(formula.terminalEvent, ~ . - 1), data.terminal)
}

# Terminal 2 Model Matrix
if((terminal2.ind == 1)) {
	if(noVarEvent[4] == 0){
		  modelmatrix4 = model.matrix(update(formula.terminalEvent2, ~ . - 1),
  		    data.terminal)
	}else{
		  modelmatrix4 = matrix(0)
	}
    terminal2 <- data.terminal[[TERMINAL2]]
} else{
    terminal2 <- data.terminal[[TERMINAL1]]*0
    modelmatrix4 = matrix(0)
}
#########################################################################
#########################################################################
# (4) Configure Parameters

### Total number of parameters
nvar = ncol(modelmatrix1) * (1-noVarEvent[1]) +
  ncol(modelmatrix2) * (1-noVarEvent[2])  +
  ncol(modelmatrix3) * event2.ind * (1-noVarEvent[3])+
  ncol(modelmatrix4) * terminal2.ind * (1-noVarEvent[4])

nbvar = c(
	ncol(modelmatrix1),
	ncol(modelmatrix2),
	ncol(modelmatrix3),
	ncol(modelmatrix4)
)

# Total number of parameters
# This will need adjustment before incorporating a second recurrent event
if(typeof==0 ){ # splines, single random effect
	np = ((2 + n.knots) * (2 + event2.ind + terminal2.ind) + nvar + 3 + 2*jointGeneral)
}else if(typeof == 2){ # weibull, single random effect
	np = (2 * (2 + event2.ind + terminal2.ind) + nvar + 3 + 2*jointGeneral)
}

#########################################################################
#########################################################################
# (5) Define GH nodes Weights

gh <- statmod::gauss.quad(GHpoints, kind="hermite")
ghNodes = gh$nodes
ghWeights = gh$weights * exp(gh$nodes^2)

############################################################
############################################################
# (6) Compute Gap Times (If Applicable)

if(gapTimes){
	tt11 <- tt11 - tt10
	tt10 <- 0 * tt10
	tt1meta0 <- tt1meta0 - tt0meta0
	tt0meta0 <- 0 * tt0meta0
}

#########################################################################
#########################################################################
# (7) Fill Parameter Vector with User-Defined Values OR Initialize Models

# Check if user entered values for hazard, input 1s if not

if(is.null(init.hazard)) init.hazard <- rep(1, np - nvar - 3 - 2*jointGeneral)

# Check if user entered values for coefficients, input 0s if not
if(is.null(init.B)) init.B <- rep(0, nvar)

# Check lengths of inputs
if(typeof == "Weibull" & length(init.hazard != 2 * (2 + event2.ind + terminal2.ind))){
	stop("init.hazard must have length 6 for weibull for three
	     event types, or length 8 for four event types.")
}else if(typeof == "Splines" & length(init.hazard != (n.knots + 2) * (2 + event2.ind + terminal2.ind))){
		stop("init.hazard must have length (n.knots + 2) * number of event types (3 or 4) for splines.")
}
if(jointGeneral & length(init.Sigma)!=3){
	stop("init.Sigma must have length 3 when jointGeneral = T.\n
	     Order should be: c(frailtySDTerminal1,  frailtySDTerminal2, frailtyCorrelation)")
}else if(!jointGeneral & length(init.Sigma)!=1){
	stop("init.Sigma must have length 1 when jointGeneral = F.")
}
if(length(init.B) != nvar){
	stop("init.B must be the same length as the number of coefficients.")
}

# If initialization desired, replace values
if(initialize){

	# ignore user-supplied initialization values if initialize == T
	init.hazard <- init.hazard*0 + 1

	init.B <- init.B*0

	# recreate time variable in original data set in case of gap times, create new formula
	if(gapTimes){
		initialization.formula <-
			paste("Surv(gapTimes, ", EVENT1, ")",
			      paste(gsub("Surv(.*)","", as.character(formula)), collapse = ""),
			      collapse = "")

		data$gapTimes <- tt11
	}else{
		initialization.formula <- formula
	}
	initialization.formula <- as.formula(initialization.formula)

	# create formula for initialization model 1 (includes terminal 1)
	initialization.formula <- terms(initialization.formula, specials = specials)
	initialization.formula1 <- drop.terms(terms(initialization.formula),
				  survival::untangle.specials(terms(initialization.formula, c("terminal2")), "terminal2", 1:10)$terms,
				  keep.response = T)
	initialization.formula1 <- formula(initialization.formula1)

	# create formula for initialization model 2 (includes terminal 2)
	initialization.formula2 <- drop.terms(terms(initialization.formula),
				  survival::untangle.specials(terms(initialization.formula, c("terminal")), "terminal", 1:10)$terms,
				  keep.response = T)
	initialization.formula2 <- formula(initialization.formula2)
	initialization.formula2 <- sub("terminal2\\(","terminal\\(",initialization.formula2)
	initialization.formula2 <- formula(paste0(initialization.formula2[2:3], collapse = "~"))


	# fit initialization model 1
	mod.joint1<-
		frailtyPenal(formula = initialization.formula1,
			 formula.terminalEvent = formula.terminalEvent,
			 jointGeneral = F,
			 data = data,
			 recurrentAG = !gapTimes,
			 hazard = "Weibull", RandDist = "LogN",
			 maxit = 100, print.times = F)

	# fit initialization model 2
	mod.joint2<-
		frailtyPenal(formula = initialization.formula2,
			 formula.terminalEvent = formula.terminalEvent2,
			 jointGeneral = F,
			 data = data,
			 recurrentAG = !gapTimes,
			 hazard = "Weibull", RandDist = "LogN",
			 maxit = 100, print.times = F)

	# Grab initialized values
		# Note: Joint model optimizes on the square root scale
		# for hazard parameters and frailty
		# variance, so we have to square to get to the original scale.
	# Recurrent Hazard
	init.hazard[1:(2+n.knots)] <- (mod.joint1$b[1:(2+n.knots)]^2 + mod.joint2$b[1:(2+n.knots)]^2)/2
		# average estimates from the two models
	# Terminal 1 Hazard
	init.hazard[(3+n.knots):(4+n.knots*2)] <- mod.joint1$b[(3+n.knots):(4+n.knots*2)]^2
	# Terminal 2 Hazard
	init.hazard[(5+n.knots*2):(6+n.knots*3)] <- mod.joint2$b[(3+n.knots):(4+n.knots*2)]^2

	# Random Effect Variance
	if(!jointGeneral){
		init.Sigma <- (abs(mod.joint1$b[5+n.knots*2]) + abs(mod.joint2$b[5+n.knots*2]))/2
			# average estimates from the two models
	}else{
		init.Sigma[1] <- abs(mod.joint1$b[5+n.knots*2])
		init.Sigma[2] <- abs(mod.joint2$b[5+n.knots*2])
		init.Sigma[3] <- 0 # rho, covariance
	}

	# Alpha
	init.Alpha1 <- mod.joint1$b[6+n.knots*2]
	init.Alpha2 <- mod.joint2$b[6+n.knots*2]

	# Coefficients
	if(noVarEvent[1] == 0){
		# average two estimates
		init.B[1:nbvar[1]] <- (mod.joint1$b[(7+n.knots*2):(6+n.knots*2+nbvar[1])] + mod.joint2$b[(7+n.knots*2):(6+n.knots*2+nbvar[1])])/2
	}
	if(noVarEvent[2] == 0){
		init.B[(1+nbvar[1]):(nbvar[1]+nbvar[2])] <- mod.joint1$b[(7+n.knots*2+nbvar[1]):(6+n.knots*2+nbvar[1]+nbvar[2])]
	}
	if(noVarEvent[4] == 0){
		init.B[(1+nbvar[1]+nbvar[2]):(nbvar[1]+nbvar[2]+nbvar[4])] <- mod.joint2$b[(7+n.knots*2+nbvar[1]):(6+n.knots*2+nbvar[1]+nbvar[4])]
	}
}

# Fill parameter vector
if(!jointGeneral){
	b <- c(sqrt(init.hazard),
	       log(init.Sigma),
	       init.Alpha1, init.Alpha2,
	       init.B)
}else{
	b <- c(sqrt(init.hazard),
	       log(init.Sigma[1:2]), # variance
	       log((init.Sigma[3]+1)/(1-init.Sigma[3])), # rho transformed using scale-logit
	       init.Alpha1, init.Alpha2,
	       init.B)
}

# save a copy of the starting value for the output
start.b <- b

if(length(b)!=np) stop("Parameter vector not the correct length.")
#debug: cat("\nmultivPenal.R:: length(b)=",length(b),
#debug:     "\nmultivPenal.R:: np = ",np,
#debug:     file='../package_tests/multiv_model_progress.dat',append=TRUE)


############################################################
############################################################
# (8) Check Dimensions of all variables for debugging
controls = c(maxit = maxit[1], # [1]
	 initialize = initialize, # [2]
	 typeof = typeof, # [3]
	 equidistant = 1, # [4]
	 irep = !crossVal, # [5] irep
	 gapTimes = gapTimes, # [6] ag0
	 nbIntervEvent = 0, # [7] nbIntervEvent
	 n.knots = n.knots, # [8]
	 event2.ind = event2.ind, # [9]
	 terminal2.ind = terminal2.ind, # [10]
	 GHpoints = GHpoints, # [11]
	 jointGeneral = as.integer(jointGeneral)) # [12] typeJoint0
if(length(controls) != 12) stop("Length of 'controls' not 12.")
#debug: cat("\nmultivPenal.R:: length(controls)=",length(controls),
#debug:     file='../package_tests/multiv_model_progress.dat',append=TRUE)

nobsEvent = c(length(event1),
	  length(terminal1),
	  length(event2))
if(length(nobsEvent)!= 3) stop("Length of 'nobsEvent' not 3.")
#debug: cat("\n\nmultivPenal.R::nobsEvent=", nobsEvent,
#debug:     file='../package_tests/multiv_model_progress.dat',append=TRUE)

if(length(kappa) != 4) stop("Length of 'kappa' not 4.")
#debug: cat("\nkappa=", kappa,file='../package_tests/multiv_model_progress.dat',append=TRUE)

### Recurrent 1
if(length(tt10) != nobsEvent[1] | length(tt10) !=  nobsEvent[1] | length(event1) !=  nobsEvent[1] | length(group1) !=  nobsEvent[1]){
	stop("Length of tt00, tt10, event1, group1 not nobsEvent[1]")
}
#debug: cat("\nmultivPenal.R:: length(tt10)=",length(tt10),
#debug:     "\nmultivPenal.R:: length(tt11)=",length(tt11),
#debug:     "\nmultivPenal.R:: length(event1)=",length(event1),
#debug:     "\nmultivPenal.R:: length(group1)=",length(group1),
#debug:     "\nmultivPenal.R:: max(group1)=",max(group1),
#debug:     "\nmultivPenal.R:: length(unique(group1))=",length(unique(group1)),
#debug:     file='../package_tests/multiv_model_progress.dat',append=TRUE)

### Recurrent 2
if(length(tt0meta0) !=  nobsEvent[3]  | length(tt1meta0) !=  nobsEvent[3]  | length(group2) !=  nobsEvent[3]  | length(event2) !=  nobsEvent[3] ){
	stop("Length of tt0meta0, tt1meta0, group2, event2 not nobsEvent[3]")
}
#debug: cat("\nmultivPenal.R:: length(tt0meta0)=",length(tt0meta0),
#debug:     "\nmultivPenal.R:: length(tt1meta0)=",length(tt1meta0),
#debug:     "\nmultivPenal.R:: length(group2)=",length(group2),
#debug:     "\nmultivPenal.R:: length(event2)=",length(event2),
#debug:     "\nmultivPenal.R:: max(group2)=",max(group2),
#debug:     "\nmultivPenal.R:: length(unique(group2))=",length(unique(group2)),
#debug:     file='../package_tests/multiv_model_progress.dat',append=TRUE)

if(length(event1) != nobsEvent[1]) stop("length(event1) != nobsEvent[1]")
if(length(event2) != nobsEvent[3]) stop("length(event2) != nobsEvent[3]")
if(length(group1) != nobsEvent[1]) stop("length(group1) != nobsEvent[1]")
if(length(group2) != nobsEvent[3]) stop("length(group2) != nobsEvent[3]")
if(length(groupdc) != nobsEvent[2]) stop("length(groupdc) != nobsEvent[2]")
if(max(group1) != nobsEvent[2]) stop("max(group1) != nobsEvent[2]")
if(max(groupdc) != nobsEvent[2]) stop("max(groupdc) != nobsEvent[2]")
if(length(tt1dc) != nobsEvent[2]) stop("length(tt1dc) != nobsEvent[2]")

### Terminal Events
if(length(terminal1) != nobsEvent[2]){
	cat("nobsEvent[2] = ", nobsEvent[2], " and length(terminal1) = ",length(terminal1))
	stop("length(terminal1) != nobsEvent[2]")
}
if(length(terminal2) != nobsEvent[2]){
	cat("nobsEvent[2] = ", nobsEvent[2], " and length(terminal2) = ",length(terminal2))
	stop("length(terminal2) != nobsEvent[2]")
}
#debug: cat("\nmultivPenal.R:: length(terminal1, terminal2, groupdc, tt1dc)=",
#debug:     length(terminal1),length(terminal2),length(groupdc),length(tt1dc),
#debug:     file='../package_tests/multiv_model_progress.dat',append=TRUE)

if(length(nbvar) != 4) stop("length(nbvar) != 4")
#debug: cat("\nmultivPenal.R:: nbvar=",nbvar,
#debug:     file='../package_tests/multiv_model_progress.dat',append=TRUE)

if(all(dim(modelmatrix1) != c(nobsEvent[1],nbvar[1]))) stop("all(dim(modelmatrix1) != c(nobsEvent[1],nbvar[1]))")
if(all(dim(modelmatrix2) != c(nobsEvent[2],nbvar[2]))) stop("all(dim(modelmatrix2) != c(nobsEvent[3],nbvar[2]))")
if(all(dim(modelmatrix3) != c(nobsEvent[3],nbvar[3]))) stop("all(dim(modelmatrix3) != c(nobsEvent[2],nbvar[3]))")
if(all(dim(modelmatrix4) != c(nobsEvent[2],nbvar[4]))) stop("all(dim(modelmatrix4) != c(nobsEvent[3],nbvar[4]))")
#debug: cat("\nmultivPenal.R:: dim(modelmatrix1, modelmatrix2, modelmatrix3, modelmatrix4)",
#debug:     dim(modelmatrix1),",", dim(modelmatrix2),",", dim(modelmatrix3),",", dim(modelmatrix4),
#debug:     file='../package_tests/multiv_model_progress.dat',append=TRUE)

if(length(noVarEvent) != 4) stop("length(noVarEvent) != 4")
#debug: cat("\nmultivPenal.R:: noVarEvent=",noVarEvent,
#debug:     file='../package_tests/multiv_model_progress.dat',append=TRUE)

if(any(is.na(modelmatrix1))|any(is.na(modelmatrix2))|any(is.na(modelmatrix3))|any(is.na(modelmatrix4))){
	stop("NA values among covariates. Reconfigure Data.")
}

######################################################################################################
# (9) Send to Fortran for optimization
    ans <- .Fortran(C_joint_multiv,
                controls = as.integer(controls),
                nobsEvent = as.integer(nobsEvent), #nobsEvent
                k0 = as.double(kappa),

                # Data Arguments
                tt00 = as.double(tt10),
                tt10 = as.double(tt11),
                tt0meta0 = as.double(tt0meta0),
                tt1meta0 = as.double(tt1meta0),
                ic0 = as.integer(event1), #ic0
                icmeta0 = as.integer(event2), #icmeta0

                groupe0 = as.integer(group1),#groupe0
                groupe0meta = as.integer(group2),#groupe0meta
                groupe0dc = as.integer(groupdc),#groupe0dc

                tt0dc0 = as.double(0*tt1dc), #tt0dc0
                tt1dc0 = as.double(tt1dc), #tt1dc0
                icdc0 = as.integer(terminal1), #icdc0
                icdc20 = as.integer(terminal2), #icdc20

                nbvar = as.integer(nbvar), #nbvar

                vax0 = as.double(modelmatrix1),
                vaxdc0 = as.double(modelmatrix2), #vaxdc0
                vaxmeta0 = as.double(modelmatrix3),
                vaxdc20 = as.double(modelmatrix4), #vaxdc20

                noVarEvent = as.integer(noVarEvent), #noVarEvent

                # Parameter Information
                np=as.integer(np),
                b=as.double(b),
                H_hessOut=as.double(matrix(0,nrow=np,ncol=np)),
                HIHOut=as.double(matrix(0,nrow=np,ncol=np)),
                resOut=as.double(0),
                LCV=as.double(rep(0,2)),
                critCV=as.integer(rep(0,6)),
                mtEvent = as.integer(rep(100,4)), #mtEvent
                mt1Event = as.integer(rep(100,4)), #mt1Event

                # Survival and Hazard Function Fits
                x1=as.double(rep(0,100)),                  #x1Out
                lam=as.double(matrix(0,nrow=100,ncol=3)),  #lamOut
                xSu1=as.double(rep(0,100)),                #xSu1
                surv=as.double(matrix(0,nrow=100,ncol=3)), #suOut
                x2=as.double(rep(0,100)),                  #x2Out
                lam2=as.double(matrix(0,nrow=100,ncol=3)), #lam2Out
                xSu2=as.double(rep(0,100)),                #xSu2
                surv2=as.double(matrix(0,nrow=100,ncol=3)),#su2Out
                x3=as.double(rep(0,100)),                  #x3Out
                lam3=as.double(matrix(0,nrow=100,ncol=3)), #lam3out
                xSu3=as.double(rep(0,100)),                #xSu3
                surv3=as.double(matrix(0,nrow=100,ncol=3)),#su3Out
                x4=as.double(rep(0,100)),                  #x4Out
    	    lam4=as.double(matrix(0,nrow=100,ncol=3)), #lam4Out
    	    xSu4=as.double(rep(0,100)),                #xSu4
    	    surv4=as.double(matrix(0,nrow=100,ncol=3)),#su4Out

                ni=as.integer(0),
                cptEvent=as.integer(rep(0,4)),
                ResMartingaleEvent=as.double(matrix(0,nrow=nobsEvent[2],ncol=3)),
                frailtyEstimates=as.double(matrix(0,nrow=nobsEvent[2],ncol=5)),

                linearpred=as.double(rep(0,nobsEvent[1])),
                linearpreddc=as.double(rep(0,nobsEvent[2])),
                linearpredM=as.double(rep(0,nobsEvent[3])),
                linearpreddc2=as.double(rep(0,nobsEvent[2])),
                ziOut1=as.double(rep(0,controls[8]+6)),
                ziOutdc=as.double(rep(0,controls[8]+6)),
                ziOutmeta=as.double(rep(0,controls[8]+6)),
                time=as.double(rep(0,controls[7]+1)),
                timedc=as.double(rep(0,controls[7]+1)),
                timeM=as.double(rep(0,controls[7]+1)),
                ghNodes = as.double(ghNodes),
                ghWeights = as.double(ghWeights),
                tolerance0 = as.double(tolerance)
    )
######################################################################################################
######################################################################################################
# (10) Format Model Summary Tables
 if(jointGeneral == F & hazard == "Weibull"){
 	f <- function(b){
 		c(b[1:6]^2,#exp(b[1:6])^2,
 		  exp(b[7]),
 		  b[(np-nvar-1):np])
 	}
 	f.prime <- function(b){
 		diag(c(2*b[1:6],#2*exp(2*b[1:6]),
 		       exp(b[7]),
 		       rep(1,nvar + 2)))
 	}
 	Parameter = c("Recurrent: Shape", "Recurrent: Scale",
 		  "Terminal1: Shape", "Terminal1: Scale",
 		  "Terminal2: Shape", "Terminal2: Scale",
 		  "Sigma",
 		  "Terminal1: Alpha", "Terminal2: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)),
 		  paste0("Terminal2: ",colnames(modelmatrix4)))
 }else if(jointGeneral == T & hazard == "Weibull"){
 	f <- function(b){
 		c(b[1:6]^2,#exp(b[1:6])^2,
 		  exp(b[7:8]),
 		  (exp(b[9]) - 1)/(exp(b[9]) + 1),
 		  b[(np-nvar-1):np])
 	}
 	f.prime <- function(b){
 		diag(c(2*b[1:6],#2*exp(2*b[1:6]),
 		       exp(b[7:8]),
 		       2*exp(2*b[9])/(1+exp(b[9]))^2,
 		       rep(1,nvar+2)))
 	}
 	Parameter = c("Recurrent: Shape", "Recurrent: Scale",
 		  "Terminal1: Shape", "Terminal1: Scale",
 		  "Terminal2: Shape", "Terminal2: Scale",
 		  "Terminal1: Sigma", "Terminal2: Sigma", "Rho",
 		  "Terminal1: Alpha", "Terminal2: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)),
 		  paste0("Terminal2: ",colnames(modelmatrix4)))
 }else if(jointGeneral == F & hazard == "Weibull"){
 	f <- function(b){
 		c(b[1:6]^2,#exp(b[1:6])^2,
 		  exp(b[7]),
 		  b[(np-nvar-1):np])
 	}
 	f.prime <- function(b){
 		diag(c(2*b[1:6],#2*exp(2*b[1:6]),
 		       exp(b[7]),
 		       rep(1,nvar + 2)))
	}
 	Parameter = c("Recurrent: Shape", "Recurrent: Scale",
 		  "Terminal1: Shape", "Terminal1: Scale",
 		  "Terminal2: Shape", "Terminal2: Scale",
 		  "Sigma",
 		  "Terminal1: Alpha", "Terminal2: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)),
 		  paste0("Terminal2: ",colnames(modelmatrix4)))
 }else if(jointGeneral == T & hazard == "Splines"){
 	f <- function(b){
 		c(b[1:6]^2,#b[1:6])^2,
 		  exp(b[7:8]),
 		  (exp(b[9]) - 1)/(exp(b[9]) + 1),
 		  b[(np-nvar-1):np])
 	}
 	f.prime <- function(b){
 		diag(c(2*b[1:6],#2*exp(2*b[1:6]),
 		       exp(b[7:8]),
 		       2*exp(2*b[9])/(1+exp(b[9]))^2,
 		       rep(1,nvar+2)))
 	}
 	Parameter = c("Recurrent: Shape", "Recurrent: Scale",
 		  "Terminal1: Shape", "Terminal1: Scale",
 		  "Terminal2: Shape", "Terminal2: Scale",
 		  "Terminal1: Sigma", "Terminal2: Sigma", "Rho",
 		  "Terminal1: Alpha", "Terminal2: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)),
 		  paste0("Terminal2: ",colnames(modelmatrix4)))
 }
 ans$varH.Raw <- matrix(ans$HIHOut, nrow = np, ncol = np)
 ans$varH.Estimate <- f.prime(ans$b) %*% ans$varH.Raw %*% f.prime(ans$b)
 ans$summary.table <- tibble(
 	Parameter = Parameter,
 	Raw = ans$b,
 	Raw.SE = sqrt(diag(ans$varH.Raw)),
 	Estimate = f(ans$b),
 	Estimate.SE = sqrt(diag(ans$varH.Estimate)),
 	LB95 = f(ans$b - 2*Raw.SE),
 	UB95 = f(ans$b + 2*Raw.SE),
 	p = 2*pnorm(q = -abs(ans$b), mean = 0, sd = Raw.SE),
 	H0 = paste(Parameter, " = ", f(rep(0, np)))
 )

######################################################################################################
######################################################################################################
# (11) Format Initialization Tables
 ans$initialization$b <- start.b

 if(initialize){
 	ans$initialization$joint1 <- mod.joint1
 	ans$initialization$joint2 <- mod.joint2

	# Extract Transformed Parameters
 	f1 <- function(b, i=3){
 		c(b[1:(length(b)-nvar+nbvar[i]-2)]^2,
 		  b[(length(b)-nvar+nbvar[i]-1):length(b)])
 	}
 	f1.prime <- function(b, i = 3){
 		diag(c(2*b[1:(length(b)-nvar+nbvar[i]-2)],
 		       rep(1,nvar-nbvar[i]+2)))
 	}
 	Parameters1 = c("Recurrent: Shape", "Recurrent: Scale",
 		  "Terminal1: Shape", "Terminal1: Scale",
 		  "Sigma", "Terminal1: Alpha",
 		  paste0("Recurrent: ",colnames(modelmatrix1)),
 		  paste0("Terminal1: ",colnames(modelmatrix2)))

 	Parameters2 = c("Recurrent: Shape", "Recurrent: Scale",
 		    "Terminal2: Shape", "Terminal2: Scale",
 		    "Sigma", "Terminal2: Alpha",
 		    paste0("Recurrent: ",colnames(modelmatrix1)),
 		    paste0("Terminal2: ",colnames(modelmatrix4)))

 	ans$initialization$varH.Estimate1 <-
 		f1.prime(ans$initialization$joint1$b) %*%
 		ans$initialization$joint1$varHtotal %*%
 		f1.prime(ans$initialization$joint1$b)
 	ans$initialization$varH.Estimate2 <-
 		f1.prime(ans$initialization$joint2$b, i=4) %*%
 		ans$initialization$joint2$varHtotal %*%
 		f1.prime(ans$initialization$joint2$b, i=4)
 	ans$initialization$summary.table1 <- tibble(
 		Parameter = Parameters1,
 		Raw = ans$initialization$joint1$b,
 		Raw.SE = sqrt(diag(ans$initialization$joint1$varHtotal)),
 		Estimate = f1(ans$initialization$joint1$b),
 		Estimate.SE = sqrt(diag(ans$initialization$varH.Estimate1)),
 		p = 2*pnorm(q = -abs(Raw), mean = 0, sd = Raw.SE)
 	)
 	ans$initialization$summary.table2 <- tibble(
 		Parameter = Parameters2,
 		Raw = ans$initialization$joint2$b,
 		Raw.SE = sqrt(diag(ans$initialization$joint2$varHtotal)),
 		Estimate = f1(ans$initialization$joint2$b, i=4),
 		Estimate.SE = sqrt(diag(ans$initialization$varH.Estimate2)),
 		LB95 = Estimate - 2*Estimate.SE,
 		UB95 = Estimate + 2*Estimate.SE,
 		p = 2*pnorm(q = -abs(Raw), mean = 0, sd = Raw.SE)
 	)

 }
return(ans)
}



