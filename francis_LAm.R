rm(list=ls())

# Funciones y Directorios ####  
library(stringr)

source('~/Documents/Rwork/Functions/Funciones/functions.R')
source('~/Documents/Rwork/Functions/Funciones/read.report.R')

dir.1<-'~/Documents/GitHub/Francis4Crabs'
dir.2<-'~/Documents/ADMwork/IFOP/2019/Lama_model/Cons_2003/norte/Lamnor2003/'
dir.3<-'~/Documents/ADMwork/IFOP/2019/Lama_model/Cons_2003/norte/A_Sens/A_Sens_NmEst/'

#unlink(dir.3,recursive=T) #borra el dir3 del Mac
#dir.create(file.path('~/Documents/ADMwork/IFOP/2019/Lama_model/Cons_2003/norte/','A_Sens/')) # crea el dir3 nuevo y vacío
#dir.create(file.path('~/Documents/ADMwork/IFOP/2019/Lama_model/Cons_2003/norte/A_Sens/','A_Sens_NmEst/')) # crea el dir3 nuevo y vacío
setwd(dir.2); file.copy(c('lamnor2003.dat', 'LAM.tpl'), dir.1)


# Corre el modelo ####
setwd(dir.3)
system('~/admb/admb LAM')
system('./LAM -ind lamnor2003.dat')


# ESTIMACIÓN TAMAÑO MUESTRAS ####

dat.file      = 'lamnor2003.dat'
data.0        <- lisread(paste(dir.3,dat.file, sep='/'));
names(data.0) <-  str_trim(names(data.0), side="right")
data.1        <- data.0
rep           <- reptoRlist('LAM.rep')                                


# Lee datos
years1   <- data.1$Ind[,1]
nyears1  <- data.1$nyrs
tallas  <- seq(10,52,1)
ntallas <- data.1$ntallas

#x  <-c(years1,rev(years1))
#x1 <-c(years1[1],years1[nyears1]+1,nyears1+1/2) #xaxp
#x2 <-c(years1[1]-1,years1[nyears1]+1) #xlim


#Proporción observada                  
pobsFm  <-rep$pobs_mflo
pobsFh  <-rep$pobs_hflo
pobsRm  <-rep$pobs_mcru
pobsRh  <-rep$pobs_hcru
#Proporción predicha
ppredFm<-rep$ppred_mflo
ppredFh<-rep$Ppred_hflo
ppredRm<-rep$ppred_mcru
ppredRh<-rep$ppred_hcru

resflm <-matrix(ncol=ntallas,nrow=nyears1)
for(i in 1:nyears1){
  for(j in 1:ntallas){
    resflm[,j]<-pobsFm[,j]-ppredFm[,j]}}

#Proporciones                                                
pFm   <- c(pobsFm,ppredFm); pFm[pFm==0]  <-NA
pFh   <- c(pobsFh,ppredFh); pFh[pFh==0]  <-NA
pRm   <- c(pobsRm,ppredRm); pRm[pRm==0]  <-NA
pRh   <- c(pobsRh,ppredRh); pRh[pRh==0]  <-NA

#arreglos                                                    
talla <- rep(gl((length(tallas)),length(years1),label=tallas),4)
años <- rep(years1,length(tallas)*4)
ind  <- c(rep("capt_obs",length(years1)*length(tallas)),
          rep("capt_est",length(years1)*length(tallas)))
pro  <- data.frame(años,talla,ind,pFm,pFh,pRm,pRh)



##----------------- Functions -----------------------------------------------------
## Ojo esto es muy especifico a como levantar los datos desde sus modelos
dw.age <- function(modelo){
  fra <- mai <- NULL
  out.jjm  <- modelo[[1]]$output[[1]] ## salidas del modelo de jurel
  in.jjm   <- modelo[[1]]$data  ## entradas del modelo Jurel
  
  tmp      <- list(effn = out.jjm[grep("EffN", names(out.jjm))],  # tamanos de muestra utilizados en las composiciones
                   phat = out.jjm[grep("phat", names(out.jjm))],  # estimacion del modelo (proporcion)
                   pobs = out.jjm[grep("pobs", names(out.jjm))],  # Observacion
                   name = c(out.jjm$Fshry_names,out.jjm$Index_names[c(1,2,4)]))  # aca van numbres de las pesquerias (ver abajO)
  
  nsam     <- nsc <- nsm <- dplyr::bind_cols(data.frame(in.jjm$Fagesample[,c(-3,-4)]), # Aca solo ordeno informacion de entrada al Jurel
                                             ## Supongo que ustedes podran adaptar desde sus modelos
                                             data.frame(fishery3=in.jjm$Flengthsample[,3]),
                                             data.frame(fishery4=in.jjm$Fagesample[,4]),
                                             data.frame(in.jjm$Iagesample[,c(1,2,4)]))
  
  for(j in 1:7) {  ## el proceso de Francis ademas de ser iterativo entre ajustes del modelos
    ## Es iterativo dentro de la funcion en si misma, por tanto, el 7 es opcional y sugerido
    ## que alcanza convergencia
    # Francis (2011) T1.8
    if (j == 3 ) {  ##  Numero de flotas que se revisara el tamano de muestras. En jurel tenemos
      ##  datos estructurados desde 4 flotas, pero use 3 porque Peru se nego a ajustar sus nm
      ## Creo que estos se los dije. Por tanto, cada j corre para cada flota a lo largo de los años
      ves    <- matrix(c(10:50)) # defino un objeto con una dimesion particular (solo para almacenar rtesultados)
    } else {
      ves    <- matrix(1:12) ## numero de edades
    }
    Obs      <- tmp$pobs[[j]][,-1] %*% ves  ## Extraigo por edad (en el caso de ustedes deberia ser por bin de talla)
    ## la proporcion observada
    ## Hago una multiplicacion matricial para definir dimensiones
    Pre      <- tmp$phat[[j]][,-1] %*% ves  ##  Lo mismo pero para la estimacion del modelo
    
    vjy      <- (tmp$phat[[j]][,-1] %*% ves^2) - (Pre^2)  ## Desde aca comienza la funcion T1.8
    tmp1     <- Obs - Pre
    nns      <- nsam[!is.na(nsam[,j]),j]
    tmp2     <- (vjy/matrix(nns))^0.5
    fra[j]   <- 1 / var(tmp1/tmp2)
    nsc[,j]  <- nsam[,j]*fra[j]
    
    # McAllister & Ianelli (1997)  ## Aca lo mismo para el otro algoritmo,
    ## que tambien es iterativo pero nunca se aplica asi en IFOP (no se porque)
    tmp1     <- apply(tmp$phat[[j]][,-1] * (1 - tmp$phat[[j]][,-1]), 1, sum)
    tmp2     <- apply((tmp$pobs[[j]][,-1] - tmp$phat[[j]][,-1])^2, 1, sum)
    ele2     <- matrix(tmp1/tmp2)
    ele1     <- 1/matrix(nns)
    mai[j]   <- t(ele1) %*% ele2
    nsm[,j]  <- nsm[,j]/nsm[,j] * mai[j]  ## McAllister y Ianelli pero ponderado
  }
  names(fra) <- names(mai) <- names(nsc) <- names(nsm) <- names(nsam) <- c(tmp$name)
  ## Aca ordeno los resultados por pesqueria (ustdedes no tendran los tmp$name,
  ## eso esta en el TPL de Jurel, asi que arreglenselas con nombres asignados directos en la clase
  ## por ejemplo, tmp$name <- c('industrial','artesanal'), en caso tengan dos flotas
  
  nsca  <- as.data.frame(t(nsc))
  names(nsca) <- c(1970:2017)
  
  nsma  <- as.data.frame(t(nsm))
  names(nsma) <- c(1970:2017)
  
  nsamm <- as.data.frame(t(nsam))
  names(nsamm) <- c(1970:2017)
  
  ## lo que resulta es un arreglo por pesqueria y año de poderador que se deben multiplicar por los nm
  
  return(list(w_francis=fra, n_macian=mai, ori=nsamm, francis=nsca, macian=nsma)) #DONE
  
  ## Traten de hacerlo y luego lo conversamos
}