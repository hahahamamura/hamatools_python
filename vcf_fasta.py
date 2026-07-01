#!/usr/bin/env python3
"""
=========================================================
VCF → FASTA HAPLOTYPES GENERATOR
=========================================================

Descrição:
Script para gerar sequências FASTA de haplótipos a partir de um VCF e um genoma
de referência. Suporta SNPs, Indels e regiões múltiplas.

Funcionalidades:
✔ BED com 3 colunas (chr start end)
✔ BED com 2 colunas (chr start → SNP 1 base)
✔ Genótipos missing (. ./.) → usa referência REAL do genoma
✔ VCF .vcf e .vcf.gz
✔ FASTA .fa/.fasta e .gz
✔ Logs visuais detalhados

Exemplos:

# Região única
python vcf_fasta.py \
-v variants.vcf.gz \
-r genome.fa \
-i chr1:1000-2000 \
-o output

# BED múltiplas regiões
python vcf_fasta.py \
-v variants.vcf \
-r genome.fa \
-b regions.bed \
-o output

Formato BED suportado:
chr1    1000    2000
chr1    1000

=========================================================
"""

import argparse
import gzip
import re
import sys
from pathlib import Path
from datetime import datetime


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")


def warn(msg):
    print(f"[WARNING] {msg}")


def error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)


def open_vcf(vcf_path):
    if not Path(vcf_path).exists():
        error(f"VCF não encontrado: {vcf_path}")
        sys.exit(1)

    if vcf_path.endswith('.gz'):
        return gzip.open(vcf_path, 'rt')
    return open(vcf_path, 'r')


def parse_fasta(fasta_path):
    if not Path(fasta_path).exists():
        error(f"FASTA não encontrado: {fasta_path}")
        sys.exit(1)

    log("Carregando genoma de referência...")

    genome = {}
    current_chr = None
    current_seq = []

    opener = gzip.open if fasta_path.endswith('.gz') else open
    mode = 'rt' if fasta_path.endswith('.gz') else 'r'

    with opener(fasta_path, mode) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chr:
                    genome[current_chr] = ''.join(current_seq).upper()
                header = line[1:].split()[0]
                current_chr = header
                current_seq = []
            else:
                current_seq.append(line)

        if current_chr:
            genome[current_chr] = ''.join(current_seq).upper()

    log(f"Genoma carregado ({len(genome)} cromossomos)\n")
    return genome


def parse_region(region_str):
    match = re.match(r'(.+):(\d+)-(\d+)$', region_str)
    if not match:
        raise ValueError(f"Formato inválido: {region_str}")

    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    return chrom, start, end


def parse_vcf_variants(vcf_path, chrom, start, end, bed_mode=False):
    """
    Extrai variantes do VCF no intervalo especificado com correções robustas:
    ✔ compatibilidade chr vs sem chr
    ✔ BED 0-based → VCF 1-based
    ✔ contagem real de variantes multialélicas
    ✔ não descarta variantes sem GT
    ✔ logs de debug opcionais
    """

    # Correção BED 0-based
    if bed_mode:
        start += 1

    variants = []
    samples = []
    variant_line_count = 0
    variant_allele_count = 0

    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith('##'):
                continue

            if line.startswith('#CHROM'):
                fields = line.strip().split('\t')
                samples = fields[9:]
                continue

            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            vcf_chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt_field = fields[4]

            # compatibilidade chr vs sem chr
            if vcf_chrom != chrom and vcf_chrom.replace("chr","") != chrom.replace("chr",""):
                continue

            if pos < start or pos > end:
                continue

            alts = alt_field.split(',')

            # captura genótipos se existirem
            genotypes = []
            if len(fields) >= 10:
                format_field = fields[8]
                if 'GT' in format_field:
                    gt_idx = format_field.split(':').index('GT')
                else:
                    gt_idx = 0

                for sample_data in fields[9:]:
                    gt = sample_data.split(':')[gt_idx]
                    genotypes.append(gt)

            variants.append({
                'pos': pos,
                'ref': ref,
                'alts': alts,
                'genotypes': genotypes
            })

            variant_line_count += 1
            variant_allele_count += len(alts)

    # logs corretos
    log(f"Amostras detectadas: {len(samples)}")
    log(f"Linhas de variantes no intervalo: {variant_line_count}")
    log(f"Variantes encontradas: {variant_allele_count}")
    

    # debug opcional
    if variants:
        for v in variants[:5]:
            log(f"Variantes: REF={v['ref']} ALT={','.join(v['alts'])}")

    return variants, samples



def get_allele(genome, chrom, pos, allele_idx, ref, alts):
    """
    Retorna a sequência do alelo para o índice GT fornecido.
    Retorna None para alelos '*' (spanning deletion).
    """
    if allele_idx == '.' or allele_idx == '':
        # genótipo missing → usa referência real do genoma
        return genome[chrom][pos-1:pos-1+len(ref)]

    if allele_idx == '0':
        # alelo referência
        return ref

    try:
        idx = int(allele_idx) - 1
        if 0 <= idx < len(alts):
            allele = alts[idx]
            # '*' = spanning deletion: esta posição já foi consumida por uma
            # deleção em variante anterior neste mesmo haplótipo.
            # Não contribui com nenhuma base; o cursor avança sobre len(ref)
            # para não duplicar a região, mas nada é inserido na sequência.
            if allele == '*':
                return None
            return allele
        else:
            # índice fora dos alts → usa referência real
            return genome[chrom][pos-1:pos-1+len(ref)]
    except (ValueError, TypeError):
        pass

    return genome[chrom][pos-1:pos-1+len(ref)]


def build_haplotype_sequence(genome, chrom, start, end, variants, hap_alleles):
    ref_seq = genome[chrom][start-1:end]
    if not variants:
        return ref_seq

    sorted_vars = sorted(zip(variants, hap_alleles), key=lambda x: x[0]['pos'])

    result = []
    current_pos = start

    for variant, allele_idx in sorted_vars:
        var_pos = variant['pos']
        ref = variant['ref']
        alts = variant['alts']

        if var_pos < current_pos:
            # Posição já consumida por variante anterior (ex: deleção ou spanning *).
            # Se for um alelo '*', não avançamos o cursor pois a deleção upstream
            # já cuidou disso. Simplesmente ignoramos a linha.
            continue

        # Preenche com referência entre a posição atual e esta variante
        if var_pos > current_pos:
            result.append(genome[chrom][current_pos-1:var_pos-1])

        allele_seq = get_allele(genome, chrom, var_pos, allele_idx, ref, alts)
        
        # Adicione isso em build_haplotype_sequence logo após get_allele:
        if isinstance(allele_seq, str) and '*' in allele_seq:
            print(f"[DEBUG] * LITERAL detectado! pos={var_pos}, allele_idx={allele_idx}, alts={alts}, allele_seq='{allele_seq}'")

        if allele_seq is None:
            # Alelo '*' (spanning deletion):
            # A deleção upstream já removeu este trecho — não inserimos nada.
            # Avançamos o cursor sobre len(ref) desta linha para não duplicar
            # a posição na próxima iteração.
            current_pos = var_pos + len(ref)
        else:
            result.append(allele_seq)
            # Para deleções reais (ex: CCCA → C), avança len(ref) para pular
            # as bases deletadas da referência.
            current_pos = var_pos + len(ref)

    # Adiciona o restante da referência após a última variante
    if current_pos <= end:
        result.append(genome[chrom][current_pos-1:end])

    return ''.join(result)


def generate_fasta(vcf, genome, region, output, region_idx=None, total_regions=None):
    chrom, start, end = parse_region(region)

    # Exibe progresso no formato (atual/total)
    if region_idx is not None and total_regions is not None:
        log(f"Processando região ({region_idx}/{total_regions}): {region}")
    else:
        log(f"Processando região: {region}")

    variants, samples = parse_vcf_variants(vcf, chrom, start, end)

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)

    outfile = output / f"{chrom}_{start}_{end}.fasta"

    with open(outfile, 'w') as out:
        for i, sample in enumerate(samples):
            gts = [v['genotypes'][i] for v in variants]

            haps = []
            for gt in gts:
                if gt in ['./.', '.|.']:
                    haps.append(('.', '.'))
                else:
                    alleles = re.split(r'[|/]', gt)
                    if len(alleles) >= 2:
                        haps.append((alleles[0], alleles[1]))
                    else:
                        haps.append(('.', '.'))

            h1 = build_haplotype_sequence(
                genome, chrom, start, end,
                variants, [h[0] for h in haps]
            )
            h2 = build_haplotype_sequence(
                genome, chrom, start, end,
                variants, [h[1] for h in haps]
            )

            out.write(f">{sample}_h1 {region}\n{h1}\n")
            out.write(f">{sample}_h2 {region}\n{h2}\n")

    log(f"FASTA gerado: {outfile}\n")
    return outfile


def parse_regions_file(bed):
    if not Path(bed).exists():
        error(f"Arquivo BED não encontrado: {bed}")
        sys.exit(1)

    regions = []

    with open(bed) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            fields = re.split(r'\s+', line)

            if len(fields) >= 3:
                regions.append(f"{fields[0]}:{fields[1]}-{fields[2]}")
            elif len(fields) == 2:
                regions.append(f"{fields[0]}:{fields[1]}-{fields[1]}")
            elif ':' in line:
                regions.append(line)

    log(f"Regiões carregadas do arquivo BED: {len(regions)}\n")
    return regions


def main():
    parser = argparse.ArgumentParser(
        description="Gerador de FASTA haplotípico a partir de VCF",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-v', '--vcf', required=True, help="Arquivo VCF (.vcf ou .vcf.gz)")
    parser.add_argument('-r', '--reference', required=True, help="Genoma FASTA (.fa/.fasta/.gz)")

    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument('-i', '--interval', help="Intervalo único chr:start-end")
    g.add_argument('-b', '--bed', help="Arquivo BED com múltiplas regiões")

    parser.add_argument('-o', '--output', required=True, help="Diretório de saída")

    args = parser.parse_args()

    log("INÍCIO DO PROCESSAMENTO")

    genome = parse_fasta(args.reference)

    if args.interval:
        regions = [args.interval]
    else:
        regions = parse_regions_file(args.bed)

    generated = []
    total_regions = len(regions)

    for idx, r in enumerate(regions, 1):
        try:
            generated.append(
                generate_fasta(
                    args.vcf,
                    genome,
                    r,
                    args.output,
                    region_idx=idx,
                    total_regions=total_regions
                )
            )
        except Exception as e:
            warn(f"Falha na região {r}: {e}")



if __name__ == '__main__':
    main()