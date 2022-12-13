### LA vs Ctrl

library(data.table)
library(ropls)
library(ggplot2)
library(pheatmap)
library(MetaboSignal)
library(dplyr)


#### load all metabolism 

SampInfo <- read.table("../../../X101SC21115155-Z01-J001/meta.tsv",header = T) %>% .[c(1:4,9:12),] %>%
  mutate(Group = factor(Group, levels = c("Control","Experiment")))
rownames(SampInfo) <- gsub("_F.*","",SampInfo$Run)

meta_intensity_neg <- read.csv("../Data/meta_intensity_neg.csv")
meta_intensity_pos <- read.csv("../Data/meta_intensity_pos.csv")

meta_intensity_neg <- meta_intensity_neg[,c(1:2,10:13,15:18)]
colnames(meta_intensity_neg)[3:10] <- gsub("neg_","",colnames(meta_intensity_neg)[3:10])
meta_intensity_pos <- meta_intensity_pos[,c(1:2,10:13,15:18)]
colnames(meta_intensity_pos)[3:10] <- gsub("pos_","",colnames(meta_intensity_pos)[3:10])




#### Differential analysis

meta_intensity <- list(neg = meta_intensity_neg, pos = meta_intensity_pos)

for(x in names(meta_intensity)) {
  dt <- meta_intensity[[x]]
  dt <- as.data.frame(dt)

  dt <- dt %>%
    mutate(FC = rowMeans(.[, 7:10]) / rowMeans(.[, 3:6]),
           log2FC = log2(FC))

  dt$Pvalue <-
    lapply(seq_len(nrow(dt)), function(x) {
      t.test(dt[x, 7:10], dt[x, 3:6])$p.value
    })

  dt$FDR <- p.adjust(dt$Pvalue, method = "BH")

  rownames(dt) <- dt$ID
  sacurine.plsda <- opls(t(dt[, c(7:10, 3:6)]), SampInfo$Group, orthoI = 0)
  dt$VIP <- getVipVn(sacurine.plsda)

  meta_intensity[[x]] <- dt
}


diff_meta_intensity <- lapply(meta_intensity, function(x){
  subset(x, abs(log2FC) > log2(1.2) & Pvalue < 0.05 & VIP > 1)
})




#### Visualization

MyheatmapMS <- function(Data, title, show_rownames = FALSE,
                        border_color = "grey60", cellheight = NA) {

  anno <- as.data.frame(SampInfo[c(5:8,1:4), c("Group")])
  rownames(anno) <- colnames(Data)
  colnames(anno) <- "Group"

  pheatmap(
    Data,
    color = colorRampPalette(c("#188AD7", "white", "#DA222C"))(100),
    annotation_col = anno,
    main = title,
    # annotation_colors = anno_cor,
    show_rownames = show_rownames,
    show_colnames = F,
    cluster_cols = F,
    cluster_rows = T,
    cellwidth = 15,
    cellheight = cellheight,
    fontsize = 9,
    annotation_names_col = F,
    border_color = border_color,
    scale = "row",
    silent = FALSE
  )
}


pdf("../Figure/All_diff_meta_heatmap_all.pdf", width = 6, height = 5)
MyheatmapMS(df_all, title = "Differential metabolism (365)", border_color = NA)
dev.off()




### pathway enrichment
```{r}

df_all_up <- do.call(rbind,diff_meta_intensity) %>% subset(.,log2FC > 0)
df_all_down <- do.call(rbind,diff_meta_intensity) %>% subset(.,log2FC < 0)
df_all_UD <- list(up = df_all_up, down = df_all_down)


df_enrich_yes <- lapply(c("up", "down"), function(d) {
  meta_neg_anno <-
    fread("../../../X101SC21115471-Z02-J001/Result-X101SC21115471-Z02-J001-B1-42/2.MetAnnotation/KEGG/meta_neg_kegg_anno.xls"
    )
  meta_pos_anno <-
    fread("../../../X101SC21115471-Z02-J001/Result-X101SC21115471-Z02-J001-B1-42/2.MetAnnotation/KEGG/meta_pos_kegg_anno.xls"
    )

  meta_pos <- rbind(meta_neg_anno, meta_pos_anno)

  pathway <-
    strsplit(meta_pos$Kegg_map, split = ";") %>% unlist %>% unique
  pathway_list <-  lapply(pathway, function(p) {
    meta_pos[grep(p, meta_pos$Kegg_map), ]$ID
  })
  names(pathway_list) <- pathway


  lapply(pathway, function(p) {
    y = length(pathway_list[[p]])
    x = intersect(df_all_UD[[d]]$ID, pathway_list[[p]]) %>% length
    n = intersect(df_all_UD[[d]]$ID, unique(unlist(pathway_list))) %>% length
    N = length(unique(unlist(pathway_list)))


    Pvalue <- phyper(x - 1, y, N - y, n, lower.tail = F)  ### Hypergeometric
    ID <- do.call(rbind, strsplit(p, split = "  ")) %>% data.frame
    colnames(ID) <- c("KEGG_map", "KEGG_name")
    cbind(
      ID,
      data.frame(
        x = x,
        y = y,
        n = n,
        N = N,
        Ratio = x/y,
        EnrichDirect = ifelse((x / n) > (y / N), "YES", "NO"),
        pvalue = Pvalue,
        ID = paste(intersect(df_all_UD[[d]]$ID, pathway_list[[p]]), collapse = ", "),
        Kegg_ID = paste(meta_pos$Kegg_ID[meta_pos$ID %in% intersect(df_all_UD[[d]]$ID, pathway_list[[p]])], collapse = ", "),
        Name = paste(meta_pos$Name[meta_pos$ID %in% intersect(df_all_UD[[d]]$ID, pathway_list[[p]])], collapse = "; ")
      )
    )
  }) %>% do.call(rbind, .) -> df_enrich

  df_enrich$padjust <- p.adjust(df_enrich$pvalue, "BH")

  df_enrich_yes <-
    subset(df_enrich, EnrichDirect == "YES") %>% .[order(.$pvalue, decreasing = F), ]

  rownames(df_enrich_yes) <- 1:nrow(df_enrich_yes)
  return(df_enrich_yes)
})

names(df_enrich_yes) <- c("up", "down")


saveRDS(df_enrich_yes, file = "../Data/KEGG_enrichment.Rds")





### heatmap-amino acid
```{r fig.height=4, fig.width=6}


up_meta_ID <- subset(df_enrich_yes[["up"]], KEGG_map %in% c("map00260","map00410","map00970","map00220",
                                              "map00350","map00270","map00340")) %>% .$ID %>%
  as.character(.) %>% strsplit(., split = ", ") %>% unlist %>% unique


down_meta_ID <- subset(df_enrich_yes[["down"]], KEGG_map %in% c("map00350","map00250","map01230","map00310","map00340","map00970",
                                                "map00290","map00471","map00220","map00520","map00400","map00280")) %>%
  .$ID %>% as.character(.) %>% strsplit(., split = ", ") %>% unlist %>% unique


pdf("../Figure/AminoAcid_meta_heatmap.pdf", width = 8, height = 4)
data <- df_all[c(up_meta_ID, down_meta_ID), ]
rownames(data) <- data$Name
data <- data[,-1]


MyheatmapMS(
  data,
  title = "Amino acid metabolism",
  border_color = NA,
  show_rownames = T
)
dev.off()


```



