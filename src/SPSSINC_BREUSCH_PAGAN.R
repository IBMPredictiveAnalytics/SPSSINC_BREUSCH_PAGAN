#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 1989, 2014
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# Author: JKP, RJC
# Version: 1.3.0

# history
# 26-Nov-2012 Adapt for R2.14; remove restriction on variance model
# 07-07-2014 html help

# To be removed...
helptext="The SPSSINC BREUSCH PAGAN command requires the R Integration Plug-in
and the R car package.

SPSSINC BREUSCH PAGAN  DEPENDENT=dependent variable 
ENTER=independent variables
[VARIANCEMODEL= list-of-predictors]
[/OPTIONS [MISSING={LISTWISE**}] [EXECUTE={TRUE**}] ]
                   {FAIL      }           {FALSE }
[/SAVE [RESIDUALSDATASET=datasetname] [COEFSDATASET=datasetname]
       [PROGRAMFILE=filespec] ]

Split files and weight are not honored by this command.

SPSSINC BREUSCH PAGAN /HELP prints this information and does nothing else.

Example:
SPSSINC BREUSCH PAGAN DEPENDENT=mpg ENTER=engine weight.

Execute the R Breusch-Pagan test for heteroscedasticity.
DEPENDENT and ENTER specify the dependent and independent
variable names.  

Categorical independent variables are automatically converted
appropriately to factors.  A constant term is automatically included.

By default, the test for heteroscedasticity is against the level of the dependent
variable.  Specify VARIANCEMODEL=v1 ... vn to test against the alternate
hypothesis that the variance is a function of variables v1 ... vn.

MISSING=LISTWISE causes listwise deletion of missing values. FAIL 
stops the procedure if missing values are encountered.

EXECUTE=FALSE runs the command syntax without running the Breusch-Pagan test.  
This is mainly useful in combination with SAVE PROGRAMFILE.

/SAVE RESIDUALSDATASET causes a dataset containing the residuals to be created.
The dataset name must not already be in use.
The case number is included as cases will only be written for input cases with no
missing data.  Data filtered by IBM SPSS Statistics are not passed to R and will not
generate cases in the residuals dataset.

COEFSDATASET causes a new dataset containing the coefficients to be created.
The dataset name must not already be in use.

PROGRAMFILE causes the R code that implements the Breusch-Pagan test to be written
to the specified file. The generated program can be a useful starting point 
for additional specifications.
"

bp<-function(dep, enter, missing="listwise", residualsdataset=NULL, coefsdataset=NULL,
       variancemodel=NULL, execute=TRUE, programfile=NULL) {

    # Run R lm and ncvTest (or older ncv.test) procedures for Breusch-Pagan heteroscedasticity test

    domain<-"SPSSINC_BREUSCH_PAGAN"
    setuplocalization(domain)
    
    gtxt <- function(...) {
	    return(gettext(...,domain=domain))
    }

    gtxtf <- function(...) {
	    return(gettextf(...,domain=domain))
    }
    tryCatch(library(car), error=function(e){
        stop(gtxtf("The R %s package is required but could not be loaded.","car"),call.=FALSE)
        }
    )

    if (identical(missing,"listwise")) {missing<-na.omit} else {missing<-na.fail}
		# Variance model variables do not need to be subset of predictors
		# but must be fetched
    allvars<-c(dep,enter)
		allvars = append(allvars, setdiff(variancemodel, allvars))
		
		
    model<-paste(dep,"~",paste(enter,collapse="+"))
    
    # ncvTest in R2.14 fails to retrieve data from this scope when it inspects the lm result
		# The workaround is to put the necessary variables in the Global environment
		# and remove them afterwards.
		
		assign(".dta", spssdata.GetDataFromSPSS(allvars,missingValueToNA=TRUE,factorMode="labels"),
				envir=.GlobalEnv)
		assign(".model", as.formula(model), envir=.GlobalEnv)
    
    res <- tryCatch(
            summary(reslm <- lm(.model, data=.dta, na.action=missing)),
            error=function(e) {return(c(gtxt("ERROR:"),e))}
           )

    if (!is.null(res$message)) {print(res$message)} else {
        miss<-ifelse(identical(missing,na.omit),"na.omit","na.fail")
		
        caption = paste(paste("lm(formula = ",model,", na.action = ",miss,")",sep=""),
            sprintf(paste("\n",gtxt("Residual standard error: "),"%.5f\n",
                        gtxt("Degrees of freedom: "),"%s\n",
                        gtxt("R-Squared: "),"%.4f\n",
                        gtxt("Adjusted R-Squared: "),"%.4f",sep=""), 
                    res$sigma, res$df[[2]], res$r.squared, res$adj.r.squared))

        coeff<-coefficients(res)
        for (i in 1:length(attributes(coeff)$dimnames[[1]])){
           attributes(coeff)$dimnames[[1]][[i]]=gtxt(attributes(coeff)$dimnames[[1]][[i]])
        }
        
        for (i in 1:length(attributes(coeff)$dimnames[[2]])){
           if(attributes(coeff)$dimnames[[2]][[i]]=="Pr(>|t|)") attributes(coeff)$dimnames[[2]][[i]] = "Sig."
           attributes(coeff)$dimnames[[2]][[i]]=gtxt(attributes(coeff)$dimnames[[2]][[i]])
        }
        
        if(packageDescription("car",fields="Version")<"2.0-0"){
		computedby = "ncv.test"
		if(!is.null(variancemodel)) {
				assign('.var.formula', as.formula(paste("~",paste(variancemodel,collapse="+"))), envir=.GlobalEnv)
				bpt = ncv.test(reslm, data=.dta, na.action=missing, var.formula = .var.formula)
		}
		else {
			trycatch(bpt = ncv.test(reslm, na.action=missing), error=function(e) print(e))
		}
	}
	else{   # data parameter is not used as of R2.14
		computedby = "ncvTest"
		if(!is.null(variancemodel)) {
					assign('.var.formula', as.formula(paste("~",paste(variancemodel,collapse="+"))), envir=.GlobalEnv)
			bpt = ncvTest(reslm, data=.dta, na.action=missing, var.formula = .var.formula)
		}
		else {

		tryCatch(bpt <- ncvTest(reslm, data=.dta, na.action=missing), error=function(e) print(e))

			}
	}
	ignore=tryCatch(remove(.dta, .model, .var.formula, envir=.GlobalVar), error=function(e) {}, warning=function(w) {})
        testresult = c(bpt$ChiSquare, bpt$Df, bpt$p)
        df = data.frame(rbind(testresult), row.names=gtxt("Test Result"))

        StartProcedure(gtxt("Residual Heteroscedasticity Test"),"SPSSINC BREUSCH PAGAN")
        
        spsspivottable.Display(coeff, 
            title=gtxt("Coefficients"), templateName="SPSSINCLM",
            caption=caption,
            isSplit=FALSE)
            
		if (is.null(variancemodel)) {
				variancemodel = gtxt("fitted values")
		} else {
				variancemodel = paste(variancemodel, collapse = "+")
		}
		caption=sprintf(paste(gtxt("Variance model: "),"%s\n",
                gtxt("Computed by R %s function"),sep=""), 
                variancemodel,
								computedby)

		collabels=c(gtxt("ChiSquare"),gtxt("D.f"), gtxt("Sig."))
		spsspivottable.Display(df,
            collabels=collabels,
            title=gtxt(bpt$test), 
			templateName="SPSSINCBREUSCHPAGAN", 
			isSplit=FALSE,
            caption= caption)
       
        spsspkg.EndProcedure()

        if (!is.null(residualsdataset)){
            dict<- spssdictionary.CreateSPSSDictionary(c("caseNumber", gtxt("Case Number"), 0, "F8.0", "nominal"),
            c("lmResiduals", model, 0, "F8.2", "scale"))
            tryCatch({
                spssdictionary.SetDictionaryToSPSS(residualsdataset, dict)
                df = data.frame(res$residuals)
                spssdata.SetDataToSPSS(residualsdataset, data.frame(row.names(df), res$residuals))
                },
                error=function(e) {print(e)
                cat(gtxt("Failed to create residuals dataset. Dataset name must not already exist: "),residualsdataset)
                }
            ) 
        }
        if (!is.null(coefsdataset)){
            dict<- spssdictionary.CreateSPSSDictionary(c("term", gtxt("Variable or Factor Value"), 100, "A100", "nominal"),
            c("coefficient", gtxt("Estimated Coefficient"), 0, "F10.3", "scale"))
            tryCatch({
                spssdictionary.SetDictionaryToSPSS(coefsdataset, dict)
                spssdata.SetDataToSPSS(coefsdataset, data.frame(row.names(res$coef), res$coef[,1]))
                },
                error=function(e) {print(e)
                cat(gtxt("Failed to create coefficients dataset. Dataset name must not already exist: "),coefsdataset)
                }
            )  
        }
        spssdictionary.EndDataStep()
    }

    res <- tryCatch(rm(list=ls()),warning=function(e){return(NULL)})
    
}

StartProcedure<-function(procname, omsid){
if (as.integer(substr(spsspkg.GetSPSSVersion(),1, 2)) >= 19)
   spsspkg.StartProcedure(procname,omsid)
else
   spsspkg.StartProcedure(omsid)
}

caller<-function(dep, enter, missing="listwise", residualsdataset=NULL, coefsdataset=NULL,
       variancemodel=NULL, programfile=NULL, execute=TRUE){

    if(!is.null(programfile)){
        lines<-c("bp<-",
            attr(bp,"source"),
            paste("dep<-",dQuote(dep),sep=""),
            paste("enter<-",deparse(enter),sep=""),
            paste("missing<-",dQuote(missing),sep=""))
        func<-"bp(dep, enter, missing"
        if(!is.null(residualsdataset)){
            func<-paste(func,", residualsdataset=",dQuote(residualsdataset),sep="")
        }
        if(!is.null(coefsdataset)){
            func<-paste(func,", coefsdataset=",dQuote(coefsdataset),sep="")
        }
        if(!is.null(variancemodel)){
            func<-paste(func,", variancemodel=",deparse(variancemodel),sep="")
        }
        func<-paste(func,")",sep="")
        lines<-c(lines,func)
        f<-file(description=programfile,open="wb",encoding="UTF-8")
        writeLines(lines,con=f)
        close(f)
    }
    
    if (execute) bp(dep, enter, missing, residualsdataset, coefsdataset, variancemodel)
    
}

setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
}  

Run<-function(args){
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
                spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep", islist=FALSE),
                spsspkg.Template("ENTER", subc="",  ktype="existingvarlist", var="enter", islist=TRUE),
                spsspkg.Template("MISSING", subc="OPTIONS",ktype="str", var="missing"),
                spsspkg.Template("VARIANCEMODEL", subc="", ktype="existingvarlist", var="variancemodel", islist=TRUE),
                spsspkg.Template("RESIDUALSDATASET", subc="SAVE", ktype="literal", var="residualsdataset"),
                spsspkg.Template("COEFSDATASET", subc="SAVE", ktype="literal", var="coefsdataset"),
                spsspkg.Template("EXECUTE", subc="OPTIONS", ktype="bool", var="execute"),
                spsspkg.Template("PROGRAMFILE", subc="SAVE", ktype="literal", var="programfile")
                ))
                
    if ("HELP" %in% attr(args,"names"))
        ###writeLines(helptext)
      helper(cmdname)
    else
        res <- spsspkg.processcmd(oobj,args,"caller")
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}