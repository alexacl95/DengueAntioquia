#####################################
### PAPER: ARBOVIRUS IN ANTIOQUIA ###
#####################################

## In this code, I do all statistical analysis from the dataset given by SSA ##

# Working Directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Loading Required Packages
library(tidyverse)
library(caret)
library(MASS)
library(car)
library(robustbase)
library(survival)
library(readxl)
library(coxrobust)
library(boot)
library(RColorBrewer)
library(viridis)
library(muhaz)
library(mice)

# Loading Datasets
datos = read_excel("Consol_deng_final.xlsx", sheet = "dataset")

#datos = datos[-c(which(datos$Region == "SIN DATO")),]

#datos$`tpo_total(hosp-inicio)` = datos$`tpo_consulta(consulta-inicio)` + datos$`tpo_deterioro(hosp-consul)`

datos[sapply(datos, is.character)] = lapply(datos[sapply(datos, is.character)], as.factor)

levels(datos$Region)[levels(datos$Region) == "URABA"] = "URABÁ"
levels(datos$Region)[levels(datos$Region) == "VALLE DE ABURRA"] = "VALLE DE ABURRÁ"

# Subsets by subregions
datos_BJ = subset(datos, Region == "BAJO CAUCA")
datos_MM = subset(datos, Region == "MAGDALENA MEDIO")
datos_ND = subset(datos, Region == "NORDESTE")
datos_N = subset(datos, Region == "NORTE")
datos_OC = subset(datos, Region == "OCCIDENTE")
datos_OR = subset(datos, Region == "ORIENTE")
datos_SO = subset(datos, Region == "SUROESTE")
datos_URA = subset(datos, Region == "URABA")
datos_VA = subset(datos, Region == "VALLE DE ABURRA")
datos_NO_DAT = subset(datos, Region == "SIN DATO")

subreg = c("datos_BJ","datos_MM","datos_ND","datos_N","datos_OC",
           "datos_OR","datos_NO_DAT","datos_SO","datos_URA","datos_VA")

variab = c(names(datos))

fun_median = function(data, idx) { 
  df = data[idx] 
  return(median(df, na.rm = TRUE))
}

#### SOCIO-DEMOGRAPHIC VARIABLES

# Variable: "EDAD"
n = length(subreg)
med.stats = matrix(0, n, 3)
impres = matrix(0,n,1)

for (i in 1:length(subreg)) {
  med.stats[i,1] = median(eval(as.symbol(paste(as.name(subreg[i]))))$edad_, na.rm = T)
  bts = boot(eval(as.symbol(paste(as.name(subreg[i]))))$edad_, fun_median, R = 1000)
  med.stats[i,2:3] = round(boot.ci(bts, conf = 0.95, type = "perc")$percent[4:5],1)
  impres[i] = paste(med.stats[i,1]," (",med.stats[i,2],"-",med.stats[i,3],")",sep = "")
}

leveneTest(edad_ ~ Region, data = datos)
kruskal.test(edad_ ~ Region, data = datos)
boxplot(edad_ ~ Region, horizontal = F, data = datos)
aaa = pairwise.wilcox.test(datos$edad_, as.factor(datos$Region), p.adjust.methods = "holm")

# Variable: "Grupos Edad"
var1 = datos$`Grupos edad`
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)
frecs.perc1 = frecs.perc

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs)

# Variable: "Sexo"
var1 = datos$sexo_
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)
frecs.perc2 = frecs.perc

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs)

# Variable: "Area"
var1 = datos$area_
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)
frecs.perc3 = frecs.perc

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs)

# Variable: "Ocupación"
var1 = datos$ocupacion_
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)
frecs.perc4 = frecs.perc

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "Pertenencia étnica"
var1 = datos$per_etn_
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)
frecs.perc5 = frecs.perc

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "Grupo Especiales"
var1 = datos$gp_discapa
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)
frecs.perc6 = frecs.perc

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

##
#### CLINICAL VARIABLES

# Variable: "Tiempo Consulta"
n = length(subreg)
med.stats = matrix(0, n, 3)
impres = matrix(0,n,1)

for (i in 1:length(subreg)) {
  med.stats[i,1] = median(eval(as.symbol(paste(as.name(subreg[i]))))$`tpo_consulta(consulta-inicio)`, na.rm = T)
  bts = boot(eval(as.symbol(paste(as.name(subreg[i]))))$`tpo_consulta(consulta-inicio)`, fun_median, R = 1000, sim = "ordinary")
  med.stats[i,2:3] = round(boot.ci(bts, conf = 0.95, type = "perc")$percent[4:5],1)
  impres[i] = paste(med.stats[i,1]," (",med.stats[i,2],"-",med.stats[i,3],")",sep = "")
}

leveneTest(`tpo_consulta(consulta-inicio)` ~ Region, data = datos)
kruskal.test(`tpo_consulta(consulta-inicio)` ~ Region, data = datos)
boxplot(`tpo_consulta(consulta-inicio)` ~ Region, horizontal = F, data = datos)
aaa = pairwise.wilcox.test(datos$`tpo_consulta(consulta-inicio)`, datos$Region, p.adjust.methods = "holm")

# Variable: "Tiempo Deterioro"
n = length(subreg)
med.stats = matrix(0, n, 3)
impres = matrix(0,n,1)

for (i in 1:length(subreg)) {
  med.stats[i,1] = median(eval(as.symbol(paste(as.name(subreg[i]))))$`tpo_deterioro(hosp-consul)`, na.rm = T)
  bts = boot(eval(as.symbol(paste(as.name(subreg[i]))))$`tpo_deterioro(hosp-consul)`, fun_median, R = 1000, sim = "ordinary")
  med.stats[i,2:3] = round(boot.ci(bts, conf = 0.95, type = "perc")$percent[4:5],1)
  impres[i] = paste(med.stats[i,1]," (",med.stats[i,2],"-",med.stats[i,3],")",sep = "")
}

leveneTest(`tpo_deterioro(hosp-consul)` ~ Region, data = datos)
kruskal.test(`tpo_deterioro(hosp-consul)` ~ Region, data = datos)
boxplot(`tpo_deterioro(hosp-consul)` ~ Region, horizontal = F, data = datos)
pairwise.wilcox.test(datos$`tpo_deterioro(hosp-consul)`, datos$Region, p.adjust.methods = "holm")

# Variable: "Hospitalizado"
var1 = datos$pac_hos_
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "Fiebre"
var1 = datos$fiebre
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 2)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "Dolor Cabeza"
var1 = datos$cefalea
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "Dolor ojo"
var1 = datos$dolrretroo
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "malgia"
var1 = datos$malgias
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "artralgia"
var1 = datos$artralgia
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "erupcion cutánea"
var1 = datos$erupcionr
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "dolor abdominal"
var1 = datos$dolor_abdo
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "vomito"
var1 = datos$vomito
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "diarrea"
var1 = datos$diarrea
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "somnolencia"
var1 = datos$somnolenci
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "hipotensión"
var1 = datos$hipotensio
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "hepatomegalia"
var1 = datos$hepatomeg
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "hematoma mucosa oral"
var1 = datos$hem_mucosa
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "hipotermia"
var1 = datos$hipotermia
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "caida plaquetas"
var1 = datos$caida_plaq
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "incrementos hematocrits"
var1 = datos$aum_hemato
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

# Variable: "dengue type"
var1 = datos$nom_eve
var2 = datos$Region
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

chisq.test(frecs.abs, correct = T)

###################################################
#### COMPARISON SEVERE-DENGUE vs NORMAL-DENGUE ####
###################################################

datos_SV = subset(datos, nom_eve == "DENGUE GRAVE")
datos_NSV = subset(datos, nom_eve == "DENGUE")
tip.deng = c("datos_SV","datos_NSV")

# Variable: "Edad"
n = length(tip.deng)
med.stats = matrix(0, n, 3)
impres = matrix(0,n,1)

for (i in 1:length(tip.deng)) {
  med.stats[i,1] = median(eval(as.symbol(paste(as.name(tip.deng[i]))))$edad_, na.rm = T)
  bts = boot(eval(as.symbol(paste(as.name(tip.deng[i]))))$edad_, fun_median, R = 1000)
  med.stats[i,2:3] = round(boot.ci(bts, conf = 0.95, type = "perc")$percent[4:5],1)
  impres[i] = paste(med.stats[i,1]," (",med.stats[i,2],"-",med.stats[i,3],")",sep = "")
}

wilcox.test(edad_ ~ nom_eve, data = datos)

# Variable: "Grupo Etáreo"
var1 = datos$`Grupos edad`
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Sexo"
var1 = datos$sexo_
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Hospitalizados"
var1 = datos$pac_hos_
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Fiebre"
var1 = datos$fiebre
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Dolor Cabeza"
var1 = datos$cefalea
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Dolor Ojo"
var1 = datos$dolrretroo
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Malgia"
var1 = datos$malgias
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Artralgia"
var1 = datos$artralgia
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Erupcion cutánea"
var1 = datos$erupcionr
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "dolor abdominal"
var1 = datos$dolor_abdo
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "vomito"
var1 = datos$vomito
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Diarrea"
var1 = datos$diarrea
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Somnolencia"
var1 = datos$somnolenci
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "hipotension"
var1 = datos$hipotensio
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Hepatomegalia"
var1 = datos$hepatomeg
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Hem mucosa"
var1 = datos$hem_mucosa
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Hipotermia"
var1 = datos$hipotermia
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Caida plaquetas"
var1 = datos$caida_plaq
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Hematocritos"
var1 = datos$aum_hemato
var2 = datos$nom_eve
n1 = length(levels(var1))
n2 = length(levels(var2))

frecs.abs = table(var1, var2)

frecs.perc = round(100 * table(var1, var2) * (1/matrix(rep(t(as.vector(table(var2))), n1), n1, n2, byrow = T)), 1)

impres = matrix(0,n1,n2)
for (i in 1:n1) {
  for (j in 1:n2) {
    impres[i,j] = paste(frecs.abs[i,j]," (",frecs.perc[i,j],"%",")",sep = "")
  }
}

m1 = nrow(datos_NSV)
m2 = nrow(datos_SV)
impres2 = matrix(0,n1,1)
for (i in 1:n1) {
  impres2[i] = prop.test(c(frecs.abs[i,]),c(m1,m2), alternative = "two.sided", correct = T)$p.value
}

# Variable: "Tiempo Consulta"
n = length(tip.deng)
med.stats = matrix(0, n, 3)
impres = matrix(0,n,1)

for (i in 1:length(tip.deng)) {
  med.stats[i,1] = median(eval(as.symbol(paste(as.name(tip.deng[i]))))$`tpo_consulta(consulta-inicio)`, na.rm = T)
  bts = boot(eval(as.symbol(paste(as.name(tip.deng[i]))))$`tpo_consulta(consulta-inicio)`, fun_median, R = 1000, sim = "ordinary")
  med.stats[i,2:3] = round(boot.ci(bts, conf = 0.95, type = "perc")$percent[4:5],1)
  impres[i] = paste(med.stats[i,1]," (",med.stats[i,2],"-",med.stats[i,3],")",sep = "")
}

wilcox.test(`tpo_consulta(consulta-inicio)` ~ nom_eve, data = datos, alternative = "t")

# Variable: "Tiempo Deterioro"
n = length(tip.deng)
med.stats = matrix(0, n, 3)
impres = matrix(0,n,1)

for (i in 1:length(subreg)) {
  med.stats[i,1] = median(eval(as.symbol(paste(as.name(tip.deng[i]))))$`tpo_deterioro(hosp-consul)`, na.rm = T)
  bts = boot(eval(as.symbol(paste(as.name(tip.deng[i]))))$`tpo_deterioro(hosp-consul)`, fun_median, R = 1000, sim = "ordinary")
  med.stats[i,2:3] = round(boot.ci(bts, conf = 0.95, type = "perc")$percent[4:5],1)
  impres[i] = paste(med.stats[i,1]," (",med.stats[i,2],"-",med.stats[i,3],")",sep = "")
}

wilcox.test(`tpo_deterioro(hosp-consul)` ~ nom_eve, data = datos, alternative = "t")

#################################### OPTIONAL !!!!!
##### SOCIO-DEMOGRAPHIC FIGURES ####
####################################

frecs.perc1 = frecs.perc1[,-c(7)]
frecs.perc2 = frecs.perc2[,-c(7)]
frecs.perc3 = frecs.perc3[,-c(7)]
frecs.perc4 = frecs.perc4[,-c(7)]
frecs.perc5 = frecs.perc5[,-c(7)]
frecs.perc6 = frecs.perc6[,-c(7)]

# FIG 1
myPalette <- brewer.pal(12, "Paired")
tiff("Img1_anex.tiff", units = "in", width = 17, height = 8, res = 600)
  par(mar = c(7,4.5,3,1))
  my_bar = barplot(frecs.perc1, beside = T, ylim = c(0,80),
                   xlab = "Subregions", main = "Age Group", ylab = "Relative Frequence (%)",
                   cex.names = 1, col = myPalette[3:8], cex.axis = 1.5)
  legend("topleft", inset = 0.02, box.lty = 1, horiz = F, text.font = 1, cex = 1, col = c(myPalette[3:8]), pch = 19, lwd = 3,
         legend = c("Adolescence","Adulthood","Childhood","Early Adulthood","Early Childhood","Old Age"))
  #text(x = my_bar, y = as.vector(frecs.perc1), label = as.vector(frecs.perc1), pos = 3, cex = 0.7, col = "black")
dev.off()
# FIG 2
tiff("Img2_anex.tiff", units = "in", width = 17, height = 8, res = 600)
  par(mar = c(7,4.5,3,1))
  my_bar = barplot(frecs.perc2, beside = T, ylim = c(0,80),
                   xlab = "Subregions", main = "Sex", ylab = "Relative Frequence (%)",
                   cex.names = 1, col = myPalette[c(6,5)], cex.axis = 1.5)
  legend("topleft", inset = 0.02, box.lty = 1, horiz = F, text.font = 1, cex = 1, col = c(myPalette[c(6,5)]), pch = 19, lwd = 3,
         legend = c("Female","Male"))
  #text(x = my_bar, y = as.vector(frecs.perc2), label = as.vector(frecs.perc2), pos = 3, cex = 0.7, col = "black")
dev.off()
# FIG 3
tiff("Img3_anex.tiff", units = "in", width = 17, height = 8, res = 600)
  par(mar = c(7,4.5,3,1))
  my_bar = barplot(frecs.perc3, beside = T, ylim = c(0,110),
                   xlab = "Subregions", main = "Type of Settlement", ylab = "Relative Frequence (%)",
                   cex.names = 1, col = myPalette[c(3,5,7)], cex.axis = 1.5)
  legend("topleft", inset = 0.02, box.lty = 1, horiz = F, text.font = 1, cex = 1, col = c(myPalette[c(3,5,7)]), pch = 19, lwd = 3,
         legend = c("Municipal Capital","Populated Center","Rural - Dispersed"))
  #text(x = my_bar, y = as.vector(frecs.perc3), label = as.vector(frecs.perc3), pos = 3, cex = 0.7, col = "black")
dev.off()
# FIG 4
tiff("Img4_anex.tiff", units = "in", width = 22, height = 9, res = 500)
  par(mar = c(7,4.5,3,1))
  my_bar = barplot(frecs.perc4, beside = T, ylim = c(0,130),
                   xlab = "Subregions", main = "Type of Occupations (ISCO-08)", ylab = "Relative Frequence (%)",
                   cex.names = 1, col = myPalette[1:10], cex.axis = 1.5)
  legend("topleft", inset = 0.02, box.lty = 1, horiz = F, text.font = 1, cex = 1, col = c(myPalette[1:10]), pch = 19, lwd = 3,
         legend = c("Skilled Agricultural, Forestry\nand Fishery Workers",
                    "Managers",
                    "Armed Forces\nOccupations",
                    "Elementary Occupations",
                    "Craft and Related\nTrades Workers",
                    "Plant and Machine\nOperators, and Assemblers",
                    "Clerical Support Workers", "Professional",
                    "Technicians and\nAssociate Professionals",
                    "Service and Sales Workers"), ncol = 2)
  #text(x = my_bar, y = as.vector(frecs.perc4), label = as.vector(frecs.perc4), pos = 3, cex = 0.7, col = "black")
dev.off()
# FIG 5
frecs.perc5 = frecs.perc5[-3,]
tiff("Img5_anex.tiff", units = "in", width = 17, height = 9, res = 500)
  par(mar = c(7,4.5,3,1))
  my_bar = barplot(frecs.perc5, beside = T, ylim = c(0,20),
                   xlab = "Subregions", main = "Ethnic Minority Groups", ylab = "Relative Frequence (%)",
                   cex.names = 1, col = myPalette[1:5], cex.axis = 1.5)
  legend("topleft", inset = 0.02, box.lty = 1, horiz = F, text.font = 1, cex = 1, col = c(myPalette[1:5]), pch = 19, lwd = 3,
         legend = c("Indigenous","Afro-Colombians\nand Mulattoes","Palenquero","Raizales","ROM"))
  #text(x = my_bar, y = as.vector(frecs.perc5), label = as.vector(frecs.perc5), pos = 3, cex = 0.7, col = "black")
dev.off()
# FIG 6
frecs.perc6 = as.data.frame(matrix(c(0.2,0.4,0.4,0.2,0.2,0.2,0.5,0.2,0.2,
                       0.7,0.6,0.3,0.2,0.5,8.5,1.2,1.9,0.2,
                       1.3,0.6,0.8,1.3,0.9,1,0.7,0.5,0.2,
                       0,0.4,0.1,0.2,0.2,0,0.3,0,0.1,
                       0.8,1.1,0.9,1.5,0.5,0.8,0.5,1,0.5,
                       0,0,0.1,0,0.1,0,0.3,0.1,0.1,
                       0,0,0,0,0.2,0,0.3,0.1,0,
                       0.1,0.4,0.2,0,0.3,1.2,0.3,0.4,0.1), 8, 9, byrow = T))
names(frecs.perc6) = c("BAJO CAUCA", "MAGDALENA MEDIO", "NORDESTE", "NORTE", "OCCIDENTE", "ORIENTE", "SUROESTE", "URABÁ", "VALLE DE ABURRÁ")
frecs.perc6 = as.matrix(frecs.perc6)
tiff("Img6_anex.tiff", units = "in", width = 17, height = 9, res = 500)
  par(mar = c(7,4.5,3,1))
  my_bar = barplot(frecs.perc6, beside = T, ylim = c(0,10),
                   xlab = "Subregions", main = "Social Grouping", ylab = "Relative Frequence (%)",
                   cex.names = 1, col = myPalette[1:8], cex.axis = 1.5)
  legend("topleft", inset = 0.02, box.lty = 1, horiz = F, text.font = 1, cex = 1, col = c(myPalette[1:8]), pch = 19, lwd = 3,
         legend = c("Disabled","Displaced","Immigrants","Expectant\nMothers","Children in\nState Care", "Demobilized","Victims of\nArmed Conflict"), ncol = 2)
  #text(x = my_bar, y = as.vector(frecs.perc6), label = as.vector(frecs.perc6), pos = 3, cex = 0.7, col = "black")
dev.off()

##############################################################
##### MODELS TO COMPARISON SEVERE-DENGUE vs NORMAL-DENGUE ####
##############################################################

# Datasets by subgroup

datos.H = subset(datos, pac_hos_ == "Sí")

#####################################
# Cox Model Dengue vs Severe Dengue #
#####################################

# Estructuras Requeridas
Status = rep(1, nrow(datos.H)) # Status: Alive (Uncensored)
Deng.Type = datos.H$nom_eve
Age.Grou = datos.H$`Grupos edad`
Sex = datos.H$sexo_
Settl.Type = datos.H$area_
Minor.Group = datos.H$per_etn_
Subr.Ant = datos.H$Region
Det_Time = datos.H$`tpo_deterioro(hosp-consul)`

d = data.frame(Det_Time, Status, Settl.Type, Sex, Age.Grou, Deng.Type, Minor.Group, Subr.Ant)
d = na.omit(d)

# Estimación de la función de Supervivencia S(t)
SurvObj = Surv(Det_Time, Status, type = "right")
S_t = survfit(SurvObj~1)
plot(S_t$surv, type = "step", lwd = 2, col = "blue")
points(S_t$surv, col = "blue", pch = 19)

# Estimación de la función de tasa de riesgo lambda(t)
Lamb_Est = -diff(log(S_t$surv))
plot(Lamb_Est, type = "step", lwd = 2, col = "blue")
points(Lamb_Est, col = "blue", pch = 19)

  Lamb_Est_2 = kphaz.fit(d$Det_Time, d$Status, method = "product-limit")$haz
  plot(Lamb_Est_2, type = "step", lwd = 2, col = "black")
  points(Lamb_Est_2, col = "black", pch = 19)
  
  plot(cumsum(Lamb_Est), type = "step", lwd = 2, col = "red")
  points(cumsum(Lamb_Est), col = "red", pch = 19)
  lines(cumsum(Lamb_Est_2), type = "step", lwd = 2, col = "black")
  points(cumsum(Lamb_Est_2), col = "black", pch = 19)

# Lanzando modelo de Cox Robusto
model1_rob = coxr(SurvObj ~ Sex + Deng.Type + Settl.Type + Subr.Ant, data = d)
model1_rob

#### Corrected version of survival curves from robust Cox model
coefficients = c(0.047, -0.104, 0.120, -0.01, 0.154, -0.153, -0.192, 0.013, -0.129, -0.073, -0.164, -0.156)
Log_Lamb_Est_cox = log(Lamb_Est)
plot(Log_Lamb_Est_cox, type = "step", lwd = 2, col = "black")
points(Log_Lamb_Est_cox, col = "black", pch = 19)

plot(S_t$time, S_t$surv, type = "step", lwd = 1, col = "white", ylim = c(0,1), xlim = c(0,12),
     xlab = "Deterioration Time (days)",
     ylab = "Empirical Survival",
     main = "Survival Functions by Sex")
grid()
lines(S_t$time, S_t_Cox_Male, type = "step", lwd = 3, col = "navy")
lines(S_t$time, S_t_Cox_Female, type = "step", lwd = 3, col = "deeppink")
legend("topright", inset = 0.02,
       legend = c("MALE","FEMALE"), col = c("navy","deeppink"),
       lty = 1, lwd = 3, box.lty = 1, horiz = FALSE, text.font = 1)

S_t_Cox_Male2 = S_t$surv ^ exp(coefficients[1] * 1)
S_t_Cox_Female2 = S_t$surv ^ exp(coefficients[1] * 0)

####

# Funciones de supervivencia Cox
S_t_Cox_Female = c(1, exp(-cumsum(exp(Log_Lamb_Est_cox + (coefficients[1] * 0)))))
S_t_Cox_Male = c(1, exp(-cumsum(exp(Log_Lamb_Est_cox + (coefficients[1] * 1)))))
#S_t_Cox_No.CentrPobl = S_t$surv ^ exp(model1_rob$coefficients[3] * 0)
#S_t_Cox_CentrPobl = S_t$surv ^ exp(model1_rob$coefficients[3] * 1)

S_t_Cox_Sub.MM = c(1, exp(-cumsum(exp(Log_Lamb_Est_cox + (coefficients[5] * 1)))))
S_t_Cox_Sub.ND = c(1, exp(-cumsum(exp(Log_Lamb_Est_cox + (coefficients[6] * 1)))))
S_t_Cox_Sub.N = c(1, exp(-cumsum(exp(Log_Lamb_Est_cox + (coefficients[7] * 1)))))
S_t_Cox_Sub.U = c(1, exp(-cumsum(exp(Log_Lamb_Est_cox + (coefficients[11] * 1)))))
S_t_Cox_Sub.VA = c(1, exp(-cumsum(exp(Log_Lamb_Est_cox + (coefficients[12] * 1)))))

#S_t_Cox_Sub.MM = S_t$surv ^ exp(model1_rob$coefficients[5] * 1)
#S_t_Cox_Sub.ND = S_t$surv ^ exp(model1_rob$coefficients[6] * 1)
#S_t_Cox_Sub.N = S_t$surv ^ exp(model1_rob$coefficients[7] * 1)
#S_t_Cox_Sub.U = S_t$surv ^ exp(model1_rob$coefficients[12] * 1)
#S_t_Cox_Sub.VA = S_t$surv ^ exp(model1_rob$coefficients[13] * 1)

# Survival functions of Robust Cox Model
#tiff("Img1.tiff", units = "in", width = 6, height = 6, res = 600)
setEPS()
postscript("Img1.eps")
  plot(S_t$time, S_t$surv, type = "step", lwd = 1, col = "white", ylim = c(0,1), xlim = c(0,12),
       xlab = "Deterioration Time (days)",
       ylab = "Empirical Survival",
       main = "Survival Functions by Sex")
  grid()
  lines(S_t$time, S_t_Cox_Male, type = "step", lwd = 3, col = "navy")
  points(S_t$time, S_t_Cox_Male, pch = 19, col = "navy")
  lines(S_t$time, S_t_Cox_Female, type = "step", lwd = 3, col = "deeppink")
  points(S_t$time, S_t_Cox_Female, pch = 19, col = "deeppink")
  legend("topright", inset = 0.02,
         legend = c("MALE","FEMALE"), col = c("navy","deeppink"),
         lty = 1, lwd = 3, box.lty = 1, horiz = FALSE, text.font = 1)
dev.off()

# tiff("Img2.tiff", units = "in", width = 7, height = 7, res = 600)
#   plot(S_t$time, S_t$surv, type = "step", lwd = 1, col = "black", ylim = c(0,1), xlim = c(0,15),
#        xlab = "Deterioration Time (days)",
#        ylab = "Empirical Survival",
#        main = "Survival Functions by Type of Seetlement")
#   grid()
#   lines(S_t$time, S_t_Cox_No.CentrPobl, type = "step", lwd = 2.5, col = "navy")
#   lines(S_t$time, S_t_Cox_CentrPobl, type = "step", lwd = 2.5, col = "deepskyblue")
#   legend("topright", inset = 0.02,
#          legend = c("Not Populated Center","Populated Center"), col = c("navy","deepskyblue"),
#          lty = 1, lwd = 3, box.lty = 1, horiz = FALSE, text.font = 1)
# dev.off()

#tiff("Img3.tiff", units = "in", width = 6, height = 6, res = 600)
setEPS()
postscript("Img3.eps")
  plot(S_t$time, S_t$surv, type = "step", lwd = 1, col = "white", ylim = c(0,1), xlim = c(0,12),
       xlab = "Deterioration Time (days)",
       ylab = "Empirical Survival",
       main = "Survival Functions by Subregion")
  grid()
  lines(S_t$time, S_t_Cox_Sub.MM, type = "step", lwd = 3, col = "chartreuse3", lty = 1)
  points(S_t$time, S_t_Cox_Sub.MM, pch = 19, col = "chartreuse3")
  lines(S_t$time, S_t_Cox_Sub.N, type = "step", lwd = 3, col = "deepskyblue", lty = 1)
  points(S_t$time, S_t_Cox_Sub.N, pch = 19, col = "deepskyblue")
  lines(S_t$time, S_t_Cox_Sub.ND, type = "step", lwd = 3, col = "navy", lty = 1)
  points(S_t$time, S_t_Cox_Sub.ND, pch = 19, col = "navy")
  lines(S_t$time, S_t_Cox_Sub.U, type = "step", lwd = 3, col = "darkorchid", lty = 1)
  points(S_t$time, S_t_Cox_Sub.U, pch = 19, col = "darkorchid")
  lines(S_t$time, S_t_Cox_Sub.VA, type = "step", lwd = 3, col = "deeppink", lty = 1)
  points(S_t$time, S_t_Cox_Sub.VA, pch = 19, col = "deeppink")
  legend("topright", inset = 0.02,
         legend = c("MAGDALENA\nMEDIO","NORTE","NORDESTE","URABÁ","VALLE\nDE ABURRÁ"),
         col = c("chartreuse3","deepskyblue","navy","darkorchid","deeppink"),
         lty = c(1,1,1,1,1), lwd = 3, box.lty = 1, horiz = FALSE, text.font = 1)
dev.off()

# Logística para modelar el tipo de dengue (normal o grave)

d.a = subset(datos, nom_eve == "DENGUE GRAVE")
d.b = subset(datos, nom_eve == "DENGUE")

NS.Indic = sample(c(1:nrow(d.b)), nrow(d.a), replace = F)
d.logMod = rbind(d.a, d.b[NS.Indic,])

#null_logit_model = glm(Deng.Type ~ 1, family = "binomial"(link = "logit"), data = d1)
logit_model = glm(nom_eve ~ `Grupos edad` + ocupacion_ + `tpo_consulta(consulta-inicio)`,
                  family = "binomial"(link = "logit"), data = d.logMod)
summary(logit_model)
step_logit_model = stepAIC(logit_model, direction = "both", trace = FALSE)
summary(step_logit_model)

# McFadden's pseudo-R squared Index
1-(logLik(logit_model)/logLik(null_logit_model)) # Logit Model
1-(logLik(step_logit_model)/logLik(null_logit_model)) # Stepwise-Logit Model


