# ------------------------------------------------------------------------
# Script for SC06 ---------------------------------------------------------
# 2018 JM update  ----------------------------------------------------------
# ------------------------------------------------------------------------

rm(list=ls())

library(R.utils)
library(jjmR)
library(tidyverse)
library(ggthemes)
library(basic.jcq)
library(scales)

sourceDirectory('../assessment/srcR')

path.input = 'input_sc06/'
path.confi = 'config_sc06/'
path.resul = 'results_sc06/'
path.fig   = 'figures_sc06_ite/'
path.exe   = '../src/build/release/jjms'
.OVERLAY   = FALSE

## Basado en la conversa de hoy, puse comentarios para que sigan el flujo

## ----------------  Reading output from base-model SC05 --------------------------
modelo     <- 'mod0.0'  ## Parte de un modelo ajustado con especificos tamanos de muestra
mod.jjm    <- readJJM(modelo, path = path.confi, input = path.input, output = path.resul)
## Esta funcion la utilizamos en la SPRFMO para leer los datos, resultados y configuracines de un modelo
## determinado por 'modelo', en este caso 'mod0.0'
## En el caso de los modelos de ustedes, podrian hacer un funcion que levante los tamanos de muestra
## y los pase a la funcion dw.age (ver comentarios en secuencia)

# iteration 1
ite.1      <- dw.age(mod.jjm)  ## Hago la iteracion -- ver funcion abajo para descripcion
                               ## Esta funcion fue ad-hoc para Jurel
modelo     <- 'mod0.4'  ## Renombro el modelo
runit(modelo, path = path.confi, input = path.input,
      est = TRUE, exec = path.exe, output = path.resul, portrait = F, pdf = TRUE)
      ## Esta funcion corre el modelo de Jurel, ustedes podrian hacer algo similar
mod.jjm.1  <- readJJM(modelo, path = path.confi, input = path.input, output = path.resul)
      ## Vuelvo a leer el modelo, pero esta vez el 'mod0.4'

# iteration 2
# Desde aca hice 5 iteraciones
# Estoy ahora viendo hacer esto en un loop, una vez tengan las suyas las discutimos para
# ver como mejorar el script
ite.2      <- dw.age(mod.jjm.1)
modelo     <- 'mod0.5'
runit(modelo, path = path.confi, input = path.input,
      est = TRUE, exec = path.exe, output = path.resul, portrait = F, pdf = TRUE)
mod.jjm.2  <- readJJM(modelo, path = path.confi, input = path.input, output = path.resul)

# iteration 3
ite.3      <- dw.age(mod.jjm.2)
modelo     <- 'mod0.6'
runit(modelo, path = path.confi, input = path.input,
      est = TRUE, exec = path.exe, output = path.resul, portrait = F, pdf = TRUE)
mod.jjm.3  <- readJJM(modelo, path = path.confi, input = path.input, output = path.resul)

# iteration 4
ite.4     <- dw.age(mod.jjm.3)
modelo     <- 'mod0.7'
runit(modelo, path = path.confi, input = path.input,
      est = TRUE, exec = path.exe, output = path.resul, portrait = F, pdf = TRUE)
mod.jjm.4  <- readJJM(modelo, path = path.confi, input = path.input, output = path.resul)

# iteration 5
ite.5     <- dw.age(mod.jjm.4)
## Este ultimo es la aplicacion de la funcion sobre el 'mod0.7'

sal <- reshape::melt(t(data.frame(ite.1$w_francis, ite.2$w_francis,
                    ite.3$w_francis, ite.4$w_francis, ite.5$w_francis)))

                    #  Esto es solo para ordenar las salidas de los diferentes modelo

## La figura que les mostre con los ponderadores (que adjunto)
ggplot(data = sal, aes(x=X1, y=value)) + geom_jitter(aes(colour=X2, group=X1), size=4, width = 0.1) +
  theme_Publication() + labs(x = '', y = '% Correction')
ggsave(file=paste0(path.fig,"all_dw",".png"),
       width = 40, height = 22, units = "cm", dpi = 600)


## ----------------  Merge several --------------------------
M <- modall <- combineModels(mod.jjm, mod.jjm.1, mod.jjm.2, mod.jjm.3, mod.jjm.4)

plot_cpue(modall, #subsetby = c("Chile_AcousN","Chile_CPUE"),
          xlab = "Year", ylab = "Survey",
          slab = "Model", ShowEstErr = TRUE, logy = FALSE, eex = 0.9)
ggsave(file=paste0(path.fig,"all_cpue",".png"),
       width = 40, height = 22, units = "cm", dpi = 600)

plot_ssb(modall, xlab = "Year", ylab = "SSB (tonnes)",
         ylim = NULL, alpha = 0.1, overlay = TRUE)
ggsave(file=paste0(path.fig,"all_ssb",".png"),
       width = 40, height = 32, units = "cm", dpi = 300)



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
