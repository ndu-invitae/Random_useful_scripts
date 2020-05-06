'''Function to extract TSS locations from GTF file, modified from GTFtools'''
def get_tss_region(GTFfile,tss_bed_file,chroms=[f'chr{x}' for x in range(1,23)]+['chrX','chrY','chrM']):
    TSSbed=[]
    with open(GTFfile,'r') as fi:
        with open(tss_bed_file,'w') as fo:
            fo.write('gene_id,gene_name,chr,strand,TSS\n') 
            for line in fi:
                table = line.split('\t')
                if table[0] in chroms and table[2] == 'gene':
                    chrom  = table[0]
                    strand = table[6]
                    geneid = line.split('gene_id')[1].split('"')[1]
                    genesymbol = line.split('gene_name')[1].split('"')[1]
                    if strand == "+":
                        fo.write(f'{geneid},{genesymbol},{chrom},{strand},{int(table[3])}\n')
                    elif strand == '-':
                        fo.write(f'{geneid},{genesymbol},{chrom},{strand},{int(table[4])}\n')
