# 1. Paper Information

## Title: Karima Echihabi, Panagiota Fatourou, Kostas Zoumpatianos, Themis Palpanas, and Houda Benbrahim. Hercules Against Data Series Similarity Search. PVLDB, 15(10), 2022. 

## Abstract: 

We propose Hercules, a parallel tree-based technique for exact similarity search on massive disk-based data series collections. We present novel index construction and query answering algorithms that leverage different summarization techniques, carefully schedule costly operations, optimize memory and disk accesses, and exploit the multi-threading and SIMD capabilities of modern hardware to perform CPU-intensive calculations. 
We demonstrate the superiority and robustness of Hercules with an extensive experimental evaluation against state-of-the-art techniques, using many synthetic and real datasets, and query workloads of varying difficulty. The results show that Hercules performs up to one order of magnitude faster than the best competitor (which is not always the same).  Moreover, Hercules is the only index that outperforms the optimized scan on all scenarios, including the hard query workloads on disk-based datasets.

