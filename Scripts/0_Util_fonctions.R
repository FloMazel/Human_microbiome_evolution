
color_host_species = c("Sorex araneus"="#8da9c4",
                       "Sorex coronatus"="#13315c",
                       "Sorex minutus" = "#a8dadc",
                       "Sorex alpinus" = "#90e0ef",
                       "Sorex antinorii" = "#00b4d8",
                       
                       "Neomys fodiens" = "#9381ff",
                       "Neomys anomalus" = "#b8b8ff",
                       
                       "Crocidura russula" = "#f9b4ed",
                       "Crocidura leucodon" = "#e574bc",
                       
                       "Talpa europaea" = "#81667a",
                  
                       
                       "Apodemus flavicollis" ="#99582a",
                       "Apodemus sylvaticus" ="#2f0e07",
                       "Apodemus alpicola" ="#bb9457",
                       
                       "Microtus agrestis"="#386641",
                       "Microtus arvalis" = "#6a994e",
                       "Microtus subterraneus" ="#a7c957",
                       
                       "Myodes glareolus"="#1b4332",
                       
                       "Chionomys nivalis"  ="#8b8c89",
                       
                       "Eliomys quercinus"  = "#d68c45",
                       
                       "Mustela nivalis" = "#f87575")

    

color_host_tissue = c("Epithelium" = "#e574bc","Feces" ="#99582a")

color_intra_inter_species_comparisons  = c("Same species" = "#3c6e71","Different species" ="#343a40")



# Define plot features 
# 1 pt = 0.35mm
# Recommended max size=7 , min size =5 in pt 
P1=5;P2=6;P3=7
MyTheme=theme_classic() + 
  theme(plot.title = element_text(size=P3),
        axis.text=element_text(size=P2),
        axis.text.x =element_text(size=P2,angle=45,hjust=1),
        axis.title=element_text(size=P2),
        
        
        
        legend.text=element_text(size=P2),
        legend.title=element_text(size=P2),
        
        panel.border = element_blank(),
        
        axis.ticks = element_line(size = 1*.35),
        axis.ticks.length = unit(.5, "mm"),
        
        plot.caption = element_text(hjust = 0), 
        plot.title.position = "plot")

MyThemeTaxaBarPlot=theme_classic() + 
  theme(plot.title = element_text(size=P3),
        axis.text=element_text(size=P2),
        axis.text.x =element_text(size=P2,angle=45,hjust=1),
        axis.title=element_text(size=P2),
        
        panel.border = element_blank(),
        
        axis.ticks = element_line(size = 1*.35),
        axis.ticks.length = unit(.5, "mm"),
        
        plot.title.position = "plot")
