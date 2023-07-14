# extraChIPs Demonstration Data

This is a repository hosting the demonstration data used in the [extraChIPs vignettes](https://www.bioconductor.org/packages/release/bioc/html/extraChIPs.html) and for the [Bioc2023](https://smped.github.io/Bioc2023_extraChIPs/) package demonstration.
Simply download the entire repository as a zip file, extract and place the complete `data` directory into your root directory when manually working through the examples.

All data was prepared from the SRA breast-cancer cell-line dataset [PRJNA509779](https://github.com/smped/PRJNA509779) using the [prepareChIPs workflow](https://github.com/smped/prepareChIPs).

All bam, narrowPeak, bed and bigwig files are exactly as taken from this workflow, but subset to the region "chr10:42354900-100000000".

## Bam Files

This was performed identically for each bam file with example code given below, before moving files to this repository

```bash
samtools view -h -b data/deduplicated/SRR8315180.sorted.bam "chr10:42354900-100000000" > data/ER/SRR8315180.bam
samtools index data/ER/SRR8315180.bam
```

## Narrow Peak Files

The same process was applied for both datasets, with example code for ERa here

```r
library(rtracklayer)
library(extraChIPs)
peakFiles <- list.files("data/ER", pattern = "narrowPeak", full.names = TRUE)
peaks <- importPeaks(peakFiles, seqinfo = sq)
subset_peaks <- endoapply(peaks, subsetByOverlaps, GRanges("chr10:42354900-100000000"))
names(subset_peaks) %>% 
  lapply(
    function(x) {
      df <- subset_peaks[[x]] %>% 
        as.data.frame %>% 
        rownames_to_column("name") %>% 
        dplyr::select(
          seqnames, start, end, name, score, strand, ends_with("Value"), peak
        ) %>%
        mutate(start = start - 1, strand = ".")
      fl <- file.path("~/github/extraChIPs_demo_data/data/ER", x)
      write.table(
        df, fl, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE
      )
    }
  )
```

## BigWig Files

The same process was applied for both datasets, with example code for ERa here

```r
library(rtracklayer)
library(extraChIPs)
fl <- list.files("output/macs2/ERa", pattern = "pileup.bw", full.names = TRUE)
gr <- GRanges("chr10:42354900-100000000")
bwfl <- BigWigFileList(fl)
names(bwfl) <- c("E2", "E2DHT")
lapply(
  names(bwfl), 
  function(x) {
    cov <- import.bw(bwfl[[x]], which = gr)
    out <- file.path(
      "~/github/extraChIPs_demo_data/data/ER", paste0(x, "_cov_chr10.bw")
    )
    export.bw(cov, out)
  }
)
fl <- list.files("output/macs2/ERa", pattern = "FE.bw", full.names = TRUE)
bwfl <- BigWigFileList(fl)
names(bwfl) <- c("E2", "E2DHT")
lapply(
  names(bwfl), 
  function(x) {
    cov <- import.bw(bwfl[[x]], which = gr)
    out <- file.path(
      "~/github/extraChIPs_demo_data/data/ER", paste0(x, "_FE_chr10.bw")
    )
    export.bw(cov, out)
  }
)
```