# AroundA
AroundA: a deep learning strategy for transcriptome-wide prediction of RNA modifications

g2t.py was used to convert genomic coordinates to transcriptomic coordinates, with the main functionality implemented via gppy.


```python g2t.py \
  --input_csv ./hallt.csv \
  --gtf_path /scratch/lb4489/bioindex/gencode.v47.annotation.gtf \
  --output_folder ./g2t_split \
  --val_chroms chr6 chr17 \
  --max_workers 20```

