import pandas as pd
import gffutils
import os
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

def init_db(path):
    global db
    import gffutils
    db = gffutils.FeatureDB(path)

def process_row(row_dict):
    chrom = row_dict["chrom"]
    pos = int(row_dict["pos"])
    strand = row_dict["strand"]
    label = row_dict["label"]

    try:
        exons = db.region(region=(chrom, pos, pos), featuretype='exon')
    except Exception:
        return []

    results = []
    for exon in exons:
        if exon.strand != strand:
            continue
        transcript_id = exon.attributes.get('transcript_id', [None])[0]
        if transcript_id:
            results.append({
                "transcript_id": transcript_id,
                "genomic_position": pos,
                "label": label
            })

    return results

def map_genomic_to_transcript(input_df, gtf_file, output_folder, db_path="gtf.db", max_workers=4):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if not os.path.exists(db_path):
        gffutils.create_db(gtf_file, dbfn=db_path, force=True, keep_order=True,
                           disable_infer_genes=True, disable_infer_transcripts=True)

    gppy_input_file = os.path.join(output_folder, "gppy_input.tsv")
    gppy_output_file = os.path.join(output_folder, "gppy_output.tsv")
    merged_output_file = os.path.join(output_folder, "gppy_output_labeled.tsv")

    print(f"\n🚀 Mapping transcript positions in parallel using {max_workers} workers...")
    
    records = []
                
    with ProcessPoolExecutor(
        max_workers=max_workers, 
        initializer=init_db, 
        initargs=(db_path,)
    ) as executor:
        futures = [
            executor.submit(process_row, row._asdict())
            for row in input_df.itertuples(index=False, name='Row')
        ]

        with tqdm(total=len(futures), desc="Mapping") as pbar:
            for future in as_completed(futures):
                result = future.result()
                if result:
                    records.extend(result)
                pbar.update(1)


    gppy_input_df = pd.DataFrame(records)
    gppy_input_df[["transcript_id", "genomic_position"]].to_csv(
        gppy_input_file, sep='\t', index=False, header=False
    )

    print("\n🛠️ Running gppy g2t...")
    
    with open(gppy_output_file, 'w') as fout:
        subprocess.run([
            "gppy", "g2t",
            "-g", gtf_file,
            "-i", gppy_input_file
        ], stdout=fout)

    g2t_df = pd.read_csv(gppy_output_file, sep='\t', header=None)
    g2t_df.columns = ["transcript_id", "genomic_position", "transcript_position", "feature"]

    final_df = pd.merge(g2t_df, gppy_input_df, on=["transcript_id", "genomic_position"], how="left")
    final_df.to_csv(merged_output_file, sep='\t', index=False)

    print(f"✅ gppy input written to: {gppy_input_file}")
    print(f"✅ labeled output written to: {merged_output_file}\n")

    return gppy_input_file, merged_output_file

def main():
    parser = argparse.ArgumentParser(description="Map genomic positions to transcript positions using gffutils and gppy")
    parser.add_argument('--input_csv', type=str, required=True, help='Input CSV with columns: chrom,pos,strand,label')
    parser.add_argument('--gtf_path', type=str, required=True, help='Path to GTF file')
    parser.add_argument('--output_folder', type=str, required=True, help='Base output directory')
    parser.add_argument('--db_path', type=str, default='gtf.db', help='Optional gffutils DB path')
    parser.add_argument('--max_workers', type=int, default=4, help='Number of parallel workers')
    parser.add_argument('--val_chroms', type=str, nargs='*', default=[], help='List of chromosomes to reserve for validation')

    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    
    input_df = pd.DataFrame({
        "chrom": df["chrom"],
        "pos": df["pos"],
        "strand": df["strand"],
        "label": df["label"] if "label" in df.columns else df["group"]
    })

    if args.val_chroms:
        val_df = input_df[input_df["chrom"].isin(args.val_chroms)]
        train_df = input_df[~input_df["chrom"].isin(args.val_chroms)]

        map_genomic_to_transcript(
            input_df=train_df,
            gtf_file=args.gtf_path,
            output_folder=os.path.join(args.output_folder, "train"),
            db_path=args.db_path,
            max_workers=args.max_workers
        )

        map_genomic_to_transcript(
            input_df=val_df,
            gtf_file=args.gtf_path,
            output_folder=os.path.join(args.output_folder, "val"),
            db_path=args.db_path,
            max_workers=args.max_workers
        )
    else:
        map_genomic_to_transcript(
            input_df=input_df,
            gtf_file=args.gtf_path,
            output_folder=args.output_folder,
            db_path=args.db_path,
            max_workers=args.max_workers
        )

if __name__ == '__main__':
    main()