rm(list = ls())

# Required packages
packages <- c("sandwich", "haven", "geigen", "pak")

# Install missing packages
missing <- packages[!vapply(packages, requireNamespace,
                            logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  install.packages(missing)
}

# Load packages
invisible(lapply(packages[1:3], library, character.only = TRUE))

# Install TwoSampleIV if needed
if (!requireNamespace("TwoSampleIV", quietly = TRUE)) {
  suppressMessages(
    pak::pkg_install("stat0357/TwoSampleIV",
                     ask = FALSE)
  )
}

# Load package
library(TwoSampleIV)




Y1 <- "dem_all"
Y2 <- "presdem_intend_all"
Y3 <- "presdem_all" #voted democrat in last election

X <- "schooling"

inst <- c("csl_16","csl_g17")
exog <- c("yearat14_2","yearat14_3","yearat14_4","yearat14_5","yearat14_6","yearat14_7",
          "yearat14_8","yearat14_9","yearat14_10","yearat14_11","yearat14_12","yearat14_13",
          "yearat14_14","yearat14_15","yearat14_16","yearat14_17","yearat14_18","yearat14_19",
          "yearat14_20","yearat14_21","yearat14_22","yearat14_23","yearat14_24","yearat14_25",
          "yearat14_26","yearat14_27","yearat14_28","yearat14_29","yearat14_30","yearat14_31",
          "yearat14_32","yearat14_33","yearat14_34","yearat14_35","yearat14_36","yearat14_37",
          "yearat14_38","yearat14_39","yearat14_40","yearat14_41","yearat14_42","yearat14_43",
          "yearat14_44","yearat14_45","yearat14_46","yearat14_47","yearat14_48","yearat14_49",
          "yearat14_50","yearat14_51","yearat14_52","yearat14_53","yearat14_54","yearat14_55",
          "yearat14_56","yearat14_57","yearat14_58","yearat14_59","yearat14_60","yearat14_61",
          "yearat14_62","yearat14_63","yearat14_64","yearat14_65","yearat14_66","yearat14_67",
          "yearat14_68","yearat14_69","yearat14_70","yearat14_71","yearat14_72","yearat14_73",
          "yearat14_74","yearat14_75","yearat14_76","yearat14_77","yearat14_78",
          "surveyyear_2","surveyyear_3","surveyyear_4","surveyyear_5","surveyyear_6",
          "state_2","state_3","state_4","state_5","state_6","state_7","state_8","state_9",
          "state_10","state_11","state_12","state_13","state_14","state_15","state_16",
          "state_17","state_18","state_19","state_20","state_21","state_22","state_23",
          "state_24","state_25","state_26","state_27","state_28","state_29","state_30",
          "state_31","state_32","state_33","state_34","state_35","state_36","state_37",
          "state_38","state_39","state_40","state_41","state_42","state_43","state_44",
          "state_45","state_46","state_47","state_48","state_49",
          "trend_2","trend_3","trend_4","trend_5","trend_6","trend_7","trend_8","trend_9",
          "trend_10","trend_11","trend_12","trend_13","trend_14","trend_15","trend_16",
          "trend_17","trend_18","trend_19","trend_20","trend_21","trend_22","trend_23",
          "trend_24","trend_25","trend_26","trend_27","trend_28","trend_29","trend_30",
          "trend_31","trend_32","trend_33","trend_34","trend_35","trend_36","trend_37",
          "trend_38","trend_39","trend_40","trend_41","trend_42","trend_43","trend_44",
          "trend_45","trend_46","trend_47","trend_48","trend_49")


dfY <- read_dta("NAES data.dta")

dfX <- read_dta("ACS data.dta")

Vhety_str = 'vcovCL(lmY,cluster=dfY$state,type="HC1",cadjust=TRUE)'

Vhetx_str = 'vcovCL(lmX,cluster=dfX$state,type="HC1",cadjust=TRUE)'



results1 <- tsiv(dfY,dfX,Y1,X,exog,inst,Vhety_str,Vhetx_str,full=FALSE,CUE=TRUE)
summary(results1)


results3 <- tsiv(dfY,dfX,Y3,X,exog,inst,Vhety_str,Vhetx_str,full=FALSE,CUE=TRUE)
summary(results3)



exog <- exog[!exog %in% c("surveyyear_2","surveyyear_5","yearat14_3","yearat14_4")]

dfY <- subset(dfY, yearat14_3==0 & yearat14_4==0)

dfX <- subset(dfX, surveyyear_2==0 & surveyyear_5==0 & yearat14_3==0 & yearat14_4==0)

results2 <- tsiv(dfY,dfX,Y2,X,exog,inst,Vhety_str,Vhetx_str,full=FALSE,CUE=TRUE)
summary(results2)

