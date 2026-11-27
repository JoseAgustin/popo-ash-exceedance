#  Programa para analizar datos en netcdf
library(ncdf4)
nc_data <- nc_open("c10M-2023.nc")
Z<- ncvar_get(nc_data,"LOAD")
dim(Z)
# Transformar el arreglo a L2
L2 <- array(0, dim = c(8760, 10000))

# Iterar sobre cada hora y reorganizar los datos en L2
for (i in 1:87600000) {
  L2[i, ] <- as.vector(Z[1, i, , ])
}

# Verificar las dimensiones de L2
dim(L2)

