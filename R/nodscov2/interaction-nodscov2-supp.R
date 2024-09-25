################################################################################
##          Support functions, dictionaries and variables
##                   to analyze interaction data
################################################################################

# Dictionaries------------------------------------------------------------------
dict_cat = c("aux nurse" = "paramedical staff", 
             "nurse" = "paramedical staff",
             "student nurse" = "paramedical staff",
             "ext physician" = "medical staff",
             "physician" = "medical staff",
             "reeducation staff" = "medical staff",
             "other" = "medical staff")

# Variables---------------------------------------------------------------------
noon_day1 = as.POSIXct(x = "2020-05-06 12:00:00")
noon_day2 = as.POSIXct(x = "2020-05-07 12:00:00")
noon_last_day = strptime("2020-08-05 12:00:00", "%Y-%m-%d %H:%M:%S")

# Functions---------------------------------------------------------------------
# Get vector of hours between two POSIXct dates 
unroll_dates = function(df) {
  data.frame(
    date_list=seq(as.POSIXct(df$DATEREMISE), as.POSIXct(df$DATEREC), "hours")
  )
}
