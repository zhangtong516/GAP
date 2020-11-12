# @Author: Zhang Tong
# @Date:   2020-11-09T16:13:56+08:00
# @Email:  zhangtong516@gmail.com
# @Filename: transcriptomic_to_genomic_coordinates.py
# @Last modified by:   Zhang Tong
# @Last modified time: 2020-11-10T11:38:24+08:00

#!/usr/bin/env python
'''
A script to transform transcript coordinates (eg: GENEA:100-200) to genomics coorditanes (eg: Chr1:150-250).
required pacagaes: GAP (https://github.com/lipan6461188/GAP.git)
@PARAMS:
    fasta file: cDNA sequence in fasta format
    gtf file: gtf file containing all transcripts
'''
##
import os, sys
import gzip
import argparse
sys.path.insert(0, '/mnt/projects/wenm/rnaStructure/ming/pipelines/icSHAPE-pipe/GAP/')
import GAP

parser = argparse.ArgumentParser(description='Convert coordinates between genomic and transcriptomic.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('mode', help='Mode of conversion. Options: genome2trans, genome2gene, gene2genome, trans2genome, gene2trans, trans2gene')
parser.add_argument('input', help='Input BED file with transcriptomic coordinates')
parser.add_argument('output', help='Output BED file with Genomic coordinates')
parser.add_argument('-c', dest='strand_column',
                    default=6, ### the format of ENCODE eCLIP bed file
                    help='The index of column which contains the strand information (+/-)')
#parser.add_argument('-f', dest='fasta',
#                    default='/mnt/projects/wenm/rnaStructure/ming/database/Genome/Ensembl/Human/GRCh38/release-100/Homo_sapiens.GRCh38.cdna.all.fa',
#                    help='Fasta file of the cDNA sequence')
parser.add_argument('-g', dest='genome_annotation',
                    default='/mnt/projects/wenm/rnaStructure/ming/database/Genome/Ensembl/Human/GRCh38/release-100/GTF/Anno.genomeCoor.bed',
                    #default='/mnt/projects/wenm/rnaStructure/ming/database/Genome/Ensembl/Human/GRCh38/release-100/Homo_sapiens.GRCh38.100.gtf',
                    help='Parsed GTF file from Ensembl-GRCh38-release100')

# parser.add_argument('-s', dest='reference_source',
#                     default='Ensembl',
#                     help='Reference source (example: Ensembl)')
parser.add_argument('-n', dest='gene_name',
                    default=False,
                    #default='/mnt/projects/wenm/rnaStructure/ming/database/Genome/Ensembl/Human/GRCh38/release-100/Homo_sapiens.GRCh38.100.gtf',
                    help='Add gene name to the output file name.')

args = parser.parse_args()

#converter = GAP.initGTF(args.gtf, genomeFile=args.fasta, source=args.reference_source)
converter = GAP.init(args.genome_annotation)
def read_bed(filein, strand_column):
    '''
    read a bed file (with or witout other information except coordinates)
    return a dictionary of coordinates
    '''
    coords = {}
    if filein.endswith("gz"):
        fin = gzip.open(filein,"rt")
    else:
        fin = open(filein, "r")

    for line in fin.readlines():
        if line[0] == '@' or line[0] == '#':
            continue
        p = line.strip().split("\t")
        chrom = p[0]
        start = int(p[1])
        end = int(p[2])
        if strand_column != 0:
            strand = p[strand_column-1]
        else:
            strand = NULL
        if len(p) > 3:
            others = ";".join( [ p[i] for i in range(3,len(p)) ] )
        else:
            others = ""
        try:
            coords[chrom].append( (chrom, start, end, strand, others))
        except KeyError:
            coords[chrom] = []
            coords[chrom].append( (chrom, start, end, strand, others))
    fin.close()
    return coords

def genome2trans(filein, fileout, strand_column):
    print("Starting converting genomic coordiates to transcriptomic coordinates...")
    coords = read_bed(filein, strand_column)

    if fileout.endswith("gz"):
        fo = gzip.open(fileout, "wt")
    else:
        fo = open(fileout,"w")

    for chrom in coords:
        print("\tConverting " + chrom + "...")
        for chr_id, start, end, strand, others in coords[chrom]:
            if chr_id.startswith("chr"):
                chr_id = chr_id[3:]
            if not strand:
                strand="+"
            try:
                covertedCoords_raw = converter.genomeCoor2transCoor(chr_id, start, end, strand)
                if len(covertedCoords_raw) == 0:
                    covertedCoords_raw = converter.genomeCoor2geneCoor(chr_id, start, end, "-")
            except:
                continue
            for region in covertedCoords_raw:
                regionID = "{}:{}-{}".format(region[0],  region[1],  region[2])
                convertedChrom = region[3]
                convertedStart = region[4]
                convertedEnd = region[5]

                try:
                    gene_name = converter.GAPer[convertedChrom]['gene_name']
                except KeyError:
                    geneInfo = converter.getGeneParser(showAttr=False)
                    gene_name = geneInfo[convertedChrom]['gene_name']
                except:
                    gene_name = NULL

                geneName_others = "{};{}".format(gene_name, others)
                fo.write("{}\t{}\t{}\t{}\t{}\n".format(convertedChrom, convertedStart, convertedEnd, regionID,
                                                       "\t".join(geneName_others.split(";"))))

    fo.close()

def genome2gene(filein, fileout, strand_column, gene_name):
    print("Starting converting genomic coordiates to transcriptomic coordinates...")
    coords = read_bed(filein, strand_column)

    if fileout.endswith("gz"):
        fo = gzip.open(fileout, "wt")
    else:
        fo = open(fileout,"w")

    for chrom in coords:
        print("\tConverting " + chrom + "...")
        for chr_id, start, end, strand, others in coords[chrom]:
            if chr_id.startswith("chr"):
                chr_id = chr_id[3:]
            if not strand:
                strand="+"
            try:
                covertedCoords_raw = converter.genomeCoor2geneCoor(chr_id, start, end, strand)
                if len(covertedCoords_raw) == 0:
                    covertedCoords_raw = converter.genomeCoor2geneCoor(chr_id, start, end, "-")
            except:
                continue
            for region in covertedCoords_raw:
                regionID = "{}:{}-{}".format(region[0],  region[1],  region[2])
                convertedChrom = region[3]
                convertedStart = region[4]
                convertedEnd = region[5]

                try:
                    geneInfo = converter.getGeneParser(showAttr=False)
                    gene_name = geneInfo[convertedChrom]['gene_name']
                except KeyError:
                    gene_name = converter.GAPer[convertedChrom]['gene_name']
                except:
                    gene_name = NULL

                geneName_others = gene_name + ";" + others
                fo.write("{}\t{}\t{}\t{}\t{}\n".format(convertedChrom, convertedStart, convertedEnd, regionID,
                                                       "\t".join(geneName_others.split(";"))))
    fo.close()

def transCoord_to_genomeCoord():
    return True
def geneCoord_to_transCoord():
    return True
def geneCoord_to_genomeCoord():
    return Ture

def main(args):
    if args.mode == "genome2trans":
        genome2trans(args.input, args.output, args.strand_column)
    if args.mode == "genome2gene":
        genome2gene(args.input, args.output, args.strand_column)

if __name__ == '__main__':
    main(args)