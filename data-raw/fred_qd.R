# ----------------
# Data transformation following the FRED-MD and FRED_QD doc.
# ----------------
trans_func = function(rawdata, tcode){

    # ==========================================================================
    # This function transforms raw data based on each series' transformation code.
    #
    # --------------------------------------------------------------------------
    # INPUT:
    #           rawdata     = raw data
    #           tcode       = transformation codes for each series
    #
    # OUTPUT:
    #           Xt          = transformed data
    # --------------------------------------------------------------------------
    # SUBFUNCTION:
    #			transxf: transforms a single series as specified
    #					 by a given transformation code.
    # ==========================================================================

    p = ncol(rawdata)

    if(p != length(tcode)) stop("The number of variables != the number of transformation codes.")

    Xt = NULL
    for(i in 1:p){

        x.temp = transxf(rawdata[, i], tcode[i])
        Xt = cbind(Xt, x.temp)

    }

    return(as.data.frame(Xt))

}


transxf = function(x, tc){

    # ==========================================================================
    # This function transforms a single series as specified
    # by a given transformation code.
    # --------------------------------------------------------------------------
    # INPUT:
    #			x 		= single series to be transformed
    #			tc  	= transformation code (1 - 7)
    # OUTPUT:
    #			xt  	= transformed series
    # ==========================================================================

    n = length(x)
    e = 1e-6
    # Allocate the output variable
    xt = rep(NA, n)

    switch(tc,
           {
               # CASE 1: Level (i.e. no transformation): x(t)
               xt = x
           },
           {
               # CASE 2: First difference: x(t)-x(t-1)
               xt[2:n] = x[2:n] - x[1:(n-1)]
           },

           {
               # CASE 3: Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
               xt[3:n] = x[3:n] - 2*x[2:(n-1)] + x[1:(n-2)]
           },
           {
               # CASE 4: Natural log: ln(x)
               if(min(x, na.rm = T) > e) xt = log(x)
           },
           {
               # CASE 5: First difference of natural log: ln(x)-ln(x-1)
               if(min(x, na.rm = T) > e){x = log(x); xt[2:n] = x[2:n] - x[1:(n-1)]}
           },
           {
               # CASE 6: Second difference of natural log: (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
               if(min(x, na.rm = T) > e){x = log(x); xt[3:n] = x[3:n] - 2*x[2:(n-1)] + x[1:(n-2)]}},
           {
               # CASE 7: First difference of percent change: (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)
               x.temp = rep(NA, n)
               x.temp[2:n] = ( x[2:n] - x[1:(n-1)] ) / x[1:(n-1)]
               xt[3:n] = x.temp[3:n] - x.temp[2:(n-1)]
           }
    )

    return(xt)

}



## code to prepare `fred_qd` dataset goes here
# A link:
# https://towardsdatascience.com/put-your-data-analysis-in-an-r-package-even-if-you-dont-publish-it-64f2bb8fd791

# Read raw data
D <- read.csv(file = "data-raw/fred_qd_data.csv")
D <- D[-nrow(D), ] # Last row is now avaiable at the time of data retriving
print(D[1:10, 1:5])
# First two rows: indicator of whether this a factor,
#   tranformation code.
# First column: date.

# Handling NAs
# For now, we simply remove all columns with missing values
# D_remove_na <- D[ , colSums(is.na(D)) == 0]

# Extract date.
date <- zoo::as.yearqtr(as.Date(D[-c(1, 2), 1], format = "%m/%d/%Y"))

# Apply transformation
tcode <- as.integer(D[2, -1])
hist(tcode) # transformation 3 and 4 are not used
D <- D[-c(1, 2), ] # We have extracted the tcode
D_trans <- trans_func(D[, -1], tcode)
colnames(D_trans) <- colnames(D)[-1]
X_trans <-subset(D_trans,
                 select = -c(UNRATE, PCECTPI, GDPCTPI, CPIAUCSL))

# Extract Inflation and unemployment rate
u <- D[, "UNRATE"] #unemployment rate (percentage point)
p_pce <- transxf(D[, "PCECTPI"], 5) * 100
p_gdp <- transxf(D[, "GDPCTPI"], 5) * 100
p_cpi <- transxf(D[, "CPIAUCSL"], 5) * 100

# Construct X matrix (y vectors are just p_pce, p_gdp, p_cpi)
QD_Phillips <- cbind(date, p_pce, p_gdp, p_cpi, unemp = u, X_trans)
QD_Phillips <- QD_Phillips[(date < 2020) & (date >= 1960), ]

usethis::use_data(QD_Phillips, overwrite = TRUE)


