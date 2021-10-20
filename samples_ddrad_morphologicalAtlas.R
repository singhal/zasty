xx = dbGetQuery(con, "select * from individuals")
x2 = dbGetQuery(con, "select * from ddrad")

xx = xx[xx$SAMPLE_ID %in% x2$SAMPLE_ID, ]
xx = cbind(xx, x2[match(xx$SAMPLE_ID, x2$SAMPLE_ID), ])
xx = xx[which(xx$GENUS == "Ctenotus"), ]

# outgroup
# alacer, ariadne, atlas, decaneurus, duricola, dux
# iapetus, impar, piankai, quatt, xenopleura
# yampiensis, zastictus, pallasotus, rhabdotus
sps = c("atlas", "ariadnae", "dux", 
        "duricola", "duricola/piankai", "grandis",
        "hanloni", "maryani", "piankai", "serventyi",
        "zastictus", "alacer", "decaneurus",
        "iapetus", "impar", "quattuordecimlineatus",
        "xenopleura", "yampiensis")
new = c("alacer", "decaneurus",
        "iapetus", "impar", "quattuordecimlineatus",
        "xenopleura", "yampiensis")
otus = c("Ctenotus_ariadnae", "Ctenotus_atlas", 
          "Ctenotus_duricola_1", "Ctenotus_duricola_2",
         "Ctenotus_dux", "Ctenotus_grandis", "Ctenotus_hanloni",
         "Ctenotus_piankai_1", "Ctenotus_piankai_2", "Ctenotus_serventyi",
         "Ctenotus_quattuordecimlineatus_2", "Ctenotus_decaneurus_1",
         "Ctenotus_decaneurus_2", "Ctenotus_iapetus", "Ctenotus_impar",
         "Ctenotus_quattuordecimlineatus_1")
newotu = c("Ctenotus_decaneurus_1",
        "Ctenotus_decaneurus_2", "Ctenotus_iapetus", "Ctenotus_impar",
        "Ctenotus_quattuordecimlineatus_1")
xx1 = xx[which(xx$SPECIES %in% sps), ]
xx2 = xx[which(xx$OTU %in% otus),]
xx3 = xx2[!xx2$SAMPLE_ID %in% xx1$SAMPLE_ID,]
xx4 = rbind(xx1, xx3)
write.csv(xx4, "~/Desktop/Ctenotus_zastictus.project.2Sept21.csv", 
          row.names = F)

# new ones only
dd= xx4[which(xx4$SPECIES %in% new | xx4$OTU %in% newotu), ]
write.csv(dd, "~/Desktop/newreads.csv", 
          row.names = F)
