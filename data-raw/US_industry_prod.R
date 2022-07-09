## code to prepare `US_industry_prod` dataset goes here
FRMD=read.csv("data-raw/FRED-MD.csv")   #independent variable
Ind_pro=read.csv("data-raw/INDPRO.csv") #dependent variable

#Data Transformation-----------------------------------------------
#(1) transform the industry index to growth rate
trans_to_rate=function(x,basis=NULL){
    if(is.null(basis)){
        basis="year"
    }
    if(basis=="year"){
        base=x[-seq(length(x)-11,length(x))]
        rate=rep(NA,12)
        rate=c(rate,diff(x,12)/base)
    } else if(basis=="month"){
        base=x[-length(x)]
        rate=NA
        rate=c(rate,diff(x,1)/base)
    }
    return(rate)
}

Ind_pro$Ind_Growth_Rate=trans_to_rate(Ind_pro$INDPRO,basis="month")
Ind_pro=Ind_pro[,-2] #drop original variable
Ind_pro=Ind_pro[-1,] #drop earliest period, whose rate is NA

#(2) change time format of FRMD (e.g., 5/1/1959->1959-05-01)
FRMD=FRMD[-1,] #discard useless rows

for(i in 1:nrow(FRMD)){
    if(nchar(FRMD[i,1])<9){
        FRMD[i,1]=paste0(0,FRMD[i,1])
    }
    FRMD[i,1]=paste(substr(FRMD[i,1],6,9),substr(FRMD[i,1],1,2),"01",sep="-")
}

#(3) merge two datasets by time index
colnames(FRMD)[1]="DATE"
US_industry_prod=merge(Ind_pro,FRMD,by="DATE",all=TRUE)

#use time index as rownames
rownames(US_industry_prod)=US_industry_prod$DATE
US_industry_prod=US_industry_prod[,-1]

#Data Cleaning (discard NAs)------------------------------------------
sum_NA=function(x){
    return(sum(is.na(x)))
}

#(1) discard periods with more than 10 NAs
NA_num_row=apply(US_industry_prod,MARGIN=1,sum_NA)
row.drop=c(seq(1,491))
US_industry_prod=US_industry_prod[-row.drop,]

#(2) discard columns with NAs
NA_num_col=apply(US_industry_prod,MARGIN=2,sum_NA)
NA_num_col=as.numeric(as.vector(NA_num_col))
US_industry_prod=US_industry_prod[,-which(NA_num_col>0)]

#Data standardization-------------------------------------------------
for(i in 1:ncol(US_industry_prod)){
    US_industry_prod[,i]=US_industry_prod[,i]/sd(US_industry_prod[,i])
}


#Save data-----------------------------------------------------
usethis::use_data(US_industry_prod, overwrite = TRUE)
