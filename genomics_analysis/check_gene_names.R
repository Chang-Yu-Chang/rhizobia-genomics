#' This script check the gene names and annotations


# checl gene names

tt$gpa %>%
    filter(str_detect(gene, "~~~")) %>%
    mutate(gene = str_remove_all(gene,"_\\d")) %>%
    mutate(gg = str_split(gene, "~~~")) %>%
    mutate(
        is_same_gene = map_lgl(gg, ~length(unique(.x)) ==1)
    ) %>%
    filter(!is_same_gene)

list_dg <- ff$per_gene_fst %>%
    # Remove unannotated genes
    filter(!str_detect(gene, "group")) %>%
    #filter(str_detect(gene, "acsA_2")) %>%
    mutate(gene = str_split(gene, "~~~")) %>%
    unnest(gene) %>%
    arrange(gene) %>%
    group_by(gene) %>%
    count %>%
    filter(n>1)



tt$gpa %>%
    filter(str_detect(gene, "group"))

sum(rowSums(tt$gpa[,-1]) == 10)
tt$gpa %>%
    filter(gene %in% tt$gpa$gene[rowSums(tt$gpa[,-1]) == 10])



# patt <- paste(list_dg$gene, collapse = "|")
# ff$per_gene_fst %>%
#     # Remove unannotated genes
#     filter(!str_detect(gene, "group")) %>%
#     filter(str_detect(gene, patt)) %>%
#     #mutate(gene = str_split(gene, "~~~")) %>%
#     #unnest(gene) %>%
#     arrange(gene) %>%
#     view
#
#
