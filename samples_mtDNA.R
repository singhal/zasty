library(RMySQL)

# the name of the database
db = 'rabosky_skinks'
# the password
pwd = 'Gorhompanpon4!'

# connect to the db
con<-dbConnect(dbDriver("MySQL"), 
               user = "raboskylab", password = pwd, 
               dbname = db, host="webapps3-db.miserver.it.umich.edu")
dbDisconnect(con)

xx = dbGetQuery(con, "select * from individuals")
x2 = dbGetQuery(con, "select * from cytbsequences")

xx = xx[xx$SAMPLE_ID %in% x2$SAMPLE_ID, ]
xx = cbind(xx, x2[match(xx$SAMPLE_ID, x2$SAMPLE_ID), ])
# most samples are identified to genus
# xx[!complete.cases(xx$GENUS), ]
xx = xx[which(xx$GENUS == "Ctenotus"), ]

# outgroup

sps = c("atlas", "ariadnae", "dux", 
        "duricola", "duricola/piankai", "grandis",
        "hanloni", "maryani", "piankai", "serventyi",
        "zastictus")
otus = c("Ctenotus_ariadnae", "Ctenotus_atlas", 
          "Ctenotus_duricola_1", "Ctenotus_duricola_2",
         "Ctenotus_dux", "Ctenotus_grandis", "Ctenotus_hanloni",
         "Ctenotus_piankai_1", "Ctenotus_piankai_2", "Ctenotus_serventyi",
         "Ctenotus_quattuordecimlineatus_2")
xx1 = xx[which(xx$SPECIES %in% sps), ]
xx2 = xx[which(xx$OTU %in% otus),]
xx3 = xx2[!xx2$SAMPLE_ID %in% xx1$SAMPLE_ID,]
xx5 = xx[xx$SAMPLE_ID %in% c("CUMV_1457414505_Ct_pant", 
        "SAMR_42153_Ct_pant", "WAMR_129335_Ct_pant"), ]
xx4 = rbind(xx1, xx3, xx5)
write.csv(xx4, "~/Desktop/Ctenotus_zastictus_mtDNA.6July21.csv", 
          row.names = F)

fileConn<-file("~/Desktop/Ctenotus_zastictus_mtDNA.cytb.fasta")
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

