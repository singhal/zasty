library(RMySQL)

# the name of the database
db = 'rabosky_skinks'
# the password
pwd = 'Gorhompanpon4!'

# connect to the db
con<-dbConnect(dbDriver("MySQL"), 
               user = "raboskylab", password = pwd, 
               dbname = db, host="webapps3-db.miserver.it.umich.edu")


xx = dbGetQuery(con, "select * from individuals")
x2 = dbGetQuery(con, "select * from cytbsequences")
dbDisconnect(con)

xx = xx[xx$SAMPLE_ID %in% x2$SAMPLE_ID, ]
xx = cbind(xx, x2[match(xx$SAMPLE_ID, x2$SAMPLE_ID), ])
# most samples are identified to genus
# xx[!complete.cases(xx$GENUS), ]
all = xx
xx = xx[which(xx$GENUS == "Ctenotus"), ]

# outgroup
sps = c("atlas", "ariadnae", "dux", 
        "duricola", "duricola/piankai", "grandis",
        "hanloni", "maryani", "piankai", "serventyi",
        "zastictus", "alacer", "decaneurus",
        "iapetus", "impar", "quattuordecimlineatus",
        "xenopleura", "yampiensis")
otus = c("Ctenotus_ariadnae", "Ctenotus_atlas", 
         "Ctenotus_duricola_1", "Ctenotus_duricola_2",
         "Ctenotus_dux", "Ctenotus_grandis", "Ctenotus_hanloni",
         "Ctenotus_piankai_1", "Ctenotus_piankai_2", "Ctenotus_serventyi",
         "Ctenotus_quattuordecimlineatus_2", "Ctenotus_decaneurus_1",
         "Ctenotus_decaneurus_2", "Ctenotus_iapetus", "Ctenotus_impar",
         "Ctenotus_quattuordecimlineatus_1")
xx1 = xx[which(xx$SPECIES %in% sps), ]
xx2 = xx[which(xx$OTU %in% otus),]
xx3 = xx2[!xx2$SAMPLE_ID %in% xx1$SAMPLE_ID,]
xx5 = all[all$SAMPLE_ID %in% c("WAMR_117238_Le_dese", 
        "SAMAR_32100_Le_dese", "CUMV_14636_Le_dese"), ]
xx4 = rbind(xx1, xx3, xx5)
write.csv(xx4, "~/Desktop/Ctenotus_zastictus_mtDNA.2Sept21.csv", 
          row.names = F)

fileConn<-file("~/Desktop/Ctenotus_zastictus_mtDNA.cytb.2Sept21.fasta")
lines = rep(NA, nrow(xx4) * 2)
for (i in 1:nrow(xx4)) {
  line1 = paste0(">", xx4[i, "SAMPLE_ID"])
  seq = xx4[i, "SEQ"]
  if (complete.cases(xx4[i, "RESEQ"])) {
    seq = xx4[i, "RESEQ"]
  }
  seq = gsub("-", "", seq)
  lines[2 * i - 1] = line1
  lines[2 * i] = seq
}
writeLines(lines, fileConn)
close(fileConn)

