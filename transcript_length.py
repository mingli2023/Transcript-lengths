import csv
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

def tss_avg_transcript_length(files, output, min_count_human, min_count_virus):
    for idx, file in enumerate(files):
        logging.info(f"Processing file: {file}")
        coordinates = {}

        with open(file, "r") as f:
            for line in f:
                line = line.strip().split("\t")
                chrom = line[0]
                start = int(line[1])
                end = int(line[2])
                ids = line[3]
                quality = int(line[4])
                strand = line[5]
                transcript_length = (end - start) + 1

                origin = start if strand == "+" else end
                coordinate_key = (chrom, origin, strand)

                if coordinate_key not in coordinates:
                    coordinates[coordinate_key] = {
                        "chr": chrom,
                        "start": origin,
                        "strand": strand,
                        "count": 1,
                        "total_transcript_length": transcript_length,
                    }
                else:
                    coordinates[coordinate_key]["count"] += 1
                    coordinates[coordinate_key]["total_transcript_length"] += transcript_length

        logging.info(f"Total unique coordinates found: {len(coordinates)}")

        for coordinate_info in coordinates.values():
            coordinate_info["average_transcript_length"] = (
                coordinate_info["total_transcript_length"] / coordinate_info["count"]
            )

        human_output = output[idx] + "-human.csv"
        virus_output = output[idx] + "-virus.csv"

        with open(human_output, "w", newline="") as csv_file_human, \
             open(virus_output, "w", newline="") as csv_file_virus:

            fieldnames = ["Chromosome", "Start", "Strand", "Count", "Average Transcript Length"]
            writer_human = csv.DictWriter(csv_file_human, fieldnames=fieldnames)
            writer_virus = csv.DictWriter(csv_file_virus, fieldnames=fieldnames)
            writer_human.writeheader()
            writer_virus.writeheader()

            for info in coordinates.values():
                if info["chr"] == "FJ616285.1" and info["count"] >= min_count_virus and info["chr"] != "JQCY02.1":
                    writer_virus.writerow({
                        "Chromosome": info["chr"],
                        "Start": info["start"],
                        "Strand": info["strand"],
                        "Count": info["count"],
                        "Average Transcript Length": info["average_transcript_length"],
                    })
                elif info["chr"] != "FJ616285.1" and info["count"] >= min_count_human and info["chr"] != "JQCY02.1":
                    writer_human.writerow({
                        "Chromosome": info["chr"],
                        "Start": info["start"],
                        "Strand": info["strand"],
                        "Count": info["count"],
                        "Average Transcript Length": info["average_transcript_length"],
                    })

        logging.info(f"Finished writing: {human_output} and {virus_output}")

if __name__ == "__main__":
    files = ["./D-NT2_Towne_WT_Flavo-dedup.bed"]
    output = ["8-dedup-tss"]
    min_count_human = 20
    min_count_virus = 10

    tss_avg_transcript_length(files, output, min_count_human, min_count_virus)
