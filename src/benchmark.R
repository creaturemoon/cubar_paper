suppressPackageStartupMessages(library(Biostrings))
library(cubar)
library(tidyverse)

yeast_cds_qc <- check_cds(yeast_cds)
yeast_cf <- count_codons(yeast_cds_qc)

yeast_exp <- yeast_exp |> 
    filter(gene_id %in% rownames(yeast_cf))

yeast_heg <- yeast_exp |> slice_max(fpkm, n = 500)
yeast_leg <- yeast_exp |> slice_min(fpkm, n = 500)

yeast_rscu <- est_rscu(yeast_cf[yeast_heg$gene_id, ])

yeast_trna_w <- est_trna_weight(yeast_trna_gcn)

yeast_opt_codons <- est_optimal_codons(
    yeast_cf[yeast_exp$gene_id, ],
    gene_score = log1p(yeast_exp$fpkm))
yeast_opt_codons <- yeast_opt_codons |>
    filter(optimal == TRUE) |>
    pull(codon)

# codon stabilization coefficient for calculating CSCg
tmp <- intersect(rownames(yeast_cf), yeast_half_life$gene_id)
yeast_half_life <- yeast_half_life |> filter(gene_id %in% tmp)
yeast_csc <- est_csc(yeast_cds_qc[yeast_half_life$gene_id], yeast_half_life)

n_seq <- c(1000, 2000, 5000, 10000, 20000, 50000, 100000,
           200000, 500000, 1000000, 2000000, 5000000, 10000000)
n_seq <- as.integer(n_seq)

mark_read_calc <- bench::press(
    n = n_seq,
    {
        bench::mark(
            min_iterations = 5, max_iterations = 10, check = FALSE,
            ENC = {
                cds <- readDNAStringSet(paste0('tmp/yeast_cds_', n, '.fasta.gz'))
                cf <- count_codons(cds)
                get_enc(cf)
            },
            CAI = {
                cds <- readDNAStringSet(paste0('tmp/yeast_cds_', n, '.fasta.gz'))
                cf <- count_codons(cds)
                get_cai(cf, rscu = yeast_rscu)
            },
            tAI = {
                cds <- readDNAStringSet(paste0('tmp/yeast_cds_', n, '.fasta.gz'))
                cf <- count_codons(cds)
                get_tai(cf, trna_w = yeast_trna_w)
            },
            Fop = {
                cds <- readDNAStringSet(paste0('tmp/yeast_cds_', n, '.fasta.gz'))
                cf <- count_codons(cds)
                get_fop(cf, op = yeast_opt_codons)
            },
            GC3s = {
                cds <- readDNAStringSet(paste0('tmp/yeast_cds_', n, '.fasta.gz'))
                cf <- count_codons(cds)
                get_gc3s(cf)
            },
            CSCg = {
                cds <- readDNAStringSet(paste0('tmp/yeast_cds_', n, '.fasta.gz'))
                cf <- count_codons(cds)
                get_cscg(cf, csc = yeast_csc)
            }
        )
    }
)
saveRDS(mark_read_calc, 'results/mark_read_calc.20240815.rds')

mark_from_seq <- bench::press(
    n = n_seq,
    {
        cds <- readDNAStringSet(paste0('tmp/yeast_cds_', n, '.fasta.gz'))
        bench::mark(
            min_iterations = 5, max_iterations = 10, check = FALSE,
            ENC = {
                cf <- count_codons(cds)
                get_enc(cf)
            },
            CAI = {
                cf <- count_codons(cds)
                get_cai(cf, rscu = yeast_rscu)
            },
            tAI = {
                cf <- count_codons(cds)
                get_tai(cf, trna_w = yeast_trna_w)
            },
            Fop = {
                cf <- count_codons(cds)
                get_fop(cf, op = yeast_opt_codons)
            },
            GC3s = {
                cf <- count_codons(cds)
                get_gc3s(cf)
            },
            CSCg = {
                cf <- count_codons(cds)
                get_cscg(cf, csc = yeast_csc)
            }
        )
    }
)
saveRDS(mark_from_seq, 'results/mark_from_seq.20240815.rds')

mark_from_cf <- bench::press(
    n = n_seq,
    {
        cds <- readDNAStringSet(paste0('tmp/yeast_cds_', n, '.fasta.gz'))
        cf <- count_codons(cds)
        bench::mark(
            min_iterations = 5, max_iterations = 10, check = FALSE,
            ENC = get_enc(cf),
            CAI = get_cai(cf, rscu = yeast_rscu),
            tAI = get_tai(cf, trna_w = yeast_trna_w),
            Fop = get_fop(cf, op = yeast_opt_codons),
            GC3s = get_gc3s(cf),
            CSCg = get_cscg(cf, csc = yeast_csc)
        )
    }
)
saveRDS(mark_from_cf, 'results/mark_from_cf.20240815.rds')
