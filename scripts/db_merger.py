def parse_fasta(file_path):
    sequences = {}
    with open(file_path) as f:
        header = ''
        sequence = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header != '':
                    sequences[sequence] = header  # Store sequence as key and header as value
                header = line[1:]
                sequence = ''
            else:
                sequence += line
        if header != '':
            sequences[sequence] = header  # Store sequence as key and header as value
    return sequences

def write_output(output_file, sequences):
    with open(output_file, 'w') as f:
        for sequence, header in sequences.items():
            f.write('>' + header + '\n')
            # Write sequence with line breaks every 60 characters
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')

def compile_databases(files):
    compiled_sequences = {}
    duplicate_sequences = {}
    for file_path in files:
        sequences = parse_fasta(file_path)
        for sequence, header in sequences.items():
            if sequence in compiled_sequences:
                if sequence not in duplicate_sequences:
                    duplicate_sequences[sequence] = [header]
                else:
                    duplicate_sequences[sequence].append(header)
            else:
                compiled_sequences[sequence] = header
    return compiled_sequences, duplicate_sequences

# List of database files
#database_files = ['./trial/AMR_CDS', './trial/db_resfinder/all.fsa']
database_files = ['./data/db/argannot/sequences',
                  './data/db/card/sequences',
                  './data/db/ecoh/sequences',
                  './data/db/ecoli_vf/sequences',
                  './data/db/megares/sequences',
                  './data/db/ncbi/sequences',
                  './data/db/plasmidfinder/sequences',
                  './data/db/resfinder/sequences',
                  './data/db/vfdb/sequences']

# Compile databases
compiled_sequences, duplicate_sequences = compile_databases(database_files)

# Write compiled sequences to a file
#write_output('compiled_AMR_database.fasta', compiled_sequences)
write_output('compiled_AMR_database.fasta', compiled_sequences)

# Write duplicate sequences to a file
# with open('abricate_duplicate_sequences.fasta', 'w') as f:
#     for sequence, headers in duplicate_sequences.items():
#         f.write(f'>Duplicate from: {", ".join(headers)}\n{sequence}\n')



