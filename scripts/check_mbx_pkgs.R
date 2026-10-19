pkgs <- c('biomformat','jsonlite','data.table','dplyr','readxl')
for (p in pkgs) { path <- system.file(package=p); cat(p,':',ifelse(nchar(path)>0,'INSTALLED','NOT INSTALLED'),'\n') }