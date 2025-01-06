import pysam

def reconstruct_sequence_with_deletions_and_quality(read):
    """Reconstruct the read sequence by adding '-' for deletions and adjusting quality."""
    sequence = []
    qualities = []
    seq_pos = 0  # Position in the original sequence
    max_quality = max(read.query_qualities) if read.query_qualities else 40

    # Iterate through the CIGAR string
    for cigar_op, length in read.cigartuples:
        if cigar_op == 0:  # Match or Mismatch (M)
            sequence.append(read.query_sequence[seq_pos:seq_pos + length])
            qualities.extend(read.query_qualities[seq_pos:seq_pos + length])
            seq_pos += length
        elif cigar_op == 2:  # Deletion (D)
            sequence.append('-' * length)
            qualities.extend([max_quality] * length)  # Use max quality for deletions
        elif cigar_op == 1:  # Insertion (I)
            seq_pos += length  # Skip the inserted bases in the query
        elif cigar_op in {4, 5}:  # Soft or Hard Clipping (S or H)
            if cigar_op == 4:  # Soft clipping affects the query sequence
                seq_pos += length

    return ''.join(sequence), qualities


def bam2fastq_subsetcontig_reconstructindel(input_bam_fn, output_fastq_fn, contig_name):
    with pysam.AlignmentFile(input_bam_fn, "rb") as bam, open(output_fastq_fn, "w") as fastq:
        # Iterate through reads mapped to the specified contig
        for read in bam.fetch(contig_name):
            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Reconstruct sequence with deletions and adjust quality
            reconstructed_sequence, reconstructed_qualities = reconstruct_sequence_with_deletions_and_quality(read)

            if reconstructed_sequence and reconstructed_qualities:
                # Convert quality scores to ASCII characters
                quality_str = ''.join(chr(q + 33) for q in reconstructed_qualities)
                fastq.write(f"@{read.query_name}\n{reconstructed_sequence}\n+\n{quality_str}\n")