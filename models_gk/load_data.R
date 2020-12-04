library(data.table)

load_data_dir = dirname(dirname(sys.frame(1)$ofile))

load_data = function() {
  btw_file = paste0(load_data_dir, "/models_gk/data/dat_btw.csv")
  wth_file = paste0(load_data_dir, "/models_gk/data/dat_wth.csv")

  dat_btw = fread(btw_file)
  dat_wth = fread(wth_file)
  
  list(dat_btw = dat_btw, dat_wth = dat_wth)
}

load_data_full = function() {
  dir_btw = paste0(load_data_dir, '/Cost2/data')
  f_btw = list.files(dir_btw, pattern=".csv", full.names=T)
  dat_btw = foreach(f = f_btw, .combine=rbind) %do% {
    info = strsplit(basename(f), "_")[[1]]
    data.table(sub=info[1], cond=info[2], date=info[3], fread(f))
  }
  dat_btw[, HandleTime := Handling + 2]
  dat_btw[, TravelTime := 20-HandleTime]
  dat_btw = dat_btw[Choice %in% c(0, 1)]
  dat_btw[, condCode := as.numeric(factor(cond, levels=c("wait", "cogTask", "phys", "pheasy")))]
  dat_btw[, half := Block > 4]
  dat_btw[, RR := sum(Offer * Choice) / sum(Choice*HandleTime + TravelTime), .(sub)]
  
  dir_wth = paste0(load_data_dir, '/Cost3/data')
  f_wth = list.files(dir_wth, pattern="main_log*.csv", full.names=T)
  dat_wth = foreach(f = f_wth, .combine=rbind) %do% {
    info = strsplit(basename(f), "_")[[1]]
    data.table(sub=info[1], cond=info[2], date=info[3], fread(f))
  }
  dat_wth[, HandleTime := Handling + 2]
  dat_wth[, TravelTime := 20-HandleTime]
  dat_wth = dat_wth[Choice %in% c(0, 1)]
  dat_wth[, cond := "wait"]
  dat_wth[Cost=="COGNITIVE", cond := "cogTask"]
  dat_wth[Cost=="GRIP", cond := "phys"]
  dat_wth[, cond := factor(cond, levels=c("wait", "cogTask", "phys"))]
  dat_wth[, condCode := as.numeric(factor(cond, levels=c("wait", "cogTask", "phys", "pheasy")))]
  dat_wth[, half := ifelse(Block > 3, "2nd", "1st")]
  dat_wth[, half := factor(half, levels=c("1st", "2nd"))]
  dat_wth[, first := ifelse(BlockType[1]==0, "cogTask", "phys"), .(sub)]
  dat_wth[, RR := sum(Offer * Choice) / sum(Choice*HandleTime + TravelTime), .(sub)]
  
  
  list(between=dat_btw, within=dat_wth)
}