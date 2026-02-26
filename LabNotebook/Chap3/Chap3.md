# Chapter 3 - Phyloseq and Rarecurve

## Purpose:
To generate the phyloseq object and rarify phyloseq object to a chosen sampling depth.

## Code:
[Phyloseq and Rarecurve](/Phyloseq_object_code.R)

## Files:
1. ms_phyloseq.Rdata (Phyloseq object)
2. ms_rare.Rdata (Rarefied sample data)

## Rarecurve:
Alpha rarefaction plot generated using R. Each line represents a distinct sample with a maximum sequencing depth of 28,263. The sampling depth was set to 6,215 as indicated by the blue dotted line.

![Rarecurve](rarecurve.png)

We selected a rarefaction depth of 6,215, where a plateau was effectively reached, and we could retain 54 samples in both PMS conditions. 
At this rarefaction depth, 3 PPMS, 1 SPMS, and 17 control samples were discarded. This selected sampling depth retained 3,424,465 (39.75%) features in 551 (96.33%) samples.