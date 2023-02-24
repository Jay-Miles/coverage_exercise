#!/usr/bin/env python

"""
Program to parse the output of sambamba, identify genes with under 30x
coverage in any exon, and produce output detailing the genes and exons
affected .

Input: text file containing sambamba output
Usage: python coverage_exercise.py <filepath>
Output: text file of unique genes; workbook of affected exon(s)

The program also uses a dump of the HGNC site to identify HGNC IDs for
genes. This file was created using the custom download tool at
https://www.genenames.org/download/custom/ using the HGNC ID, approved
symbol, previous symbols, alias symbols, and RefSeq ID columns.

Author: Jay_Miles
Date: 2023-02-24
"""


import pandas as pd
import re
import sys


def create_df(input_file):
    """ Read a sambamba output file into a pandas dataframe.

    args:
        input_file [str]: path to sambamba output text file

    returns:
        output_prefix [str]: prefix for output files
        coverage_df [df]: data converted to pandas dataframe
    """

    coverage_df = pd.read_csv(input_file, delim_whitespace=True)

    # split 'GeneSymbol;Accession' into two columns

    coverage_df[['GeneSymbol', 'Accession']] = coverage_df[
        'GeneSymbol;Accession'].str.split(';', n=1, expand=True)

    # generate the prefix for output files

    output_prefix = re.search(
        '(.*?)(?=.sambamba_output.txt$)', input_file).group()

    return output_prefix, coverage_df


def find_coverage_issues(coverage_df):
    """ Identify rows of sambamba output where x30 coverage for that
    exon is <100%, and use these to create a new dataframe.

    args:
        coverage_df [df]: pandas dataframe of sambamba output

    returns:
        issues_df [df]: subset of data for exons with 30x < 100%
    """

    issues_df = pd.DataFrame({
        'gene': [], 'accession': [], 'exon_position': [], 'percent_x30': []})

    for index, row in coverage_df.iterrows():
        if row['percentage30'] != 100:

            new_row = pd.DataFrame([[
                row['GeneSymbol'], row['Accession'],
                row['FullPosition'], row['percentage30']]],
                columns=list(issues_df))

            issues_df = pd.concat([issues_df, new_row])

    return issues_df


def get_hgncs(hgnc_dump, issues_df):
    """ Identify HGNC IDs for genes with low coverage.

    args:
        hgnc_dump [str]: path to file dump of hgnc site
        issues_df [df]: data for exons with low coverage

    returns:
        issues_df [df]: updated with gene hgnc ids
    """

    hgncs = []
    hgnc_df = pd.read_csv(hgnc_dump, sep='\t')
    hgnc_df.columns = ['hgnc', 'approved', 'previous', 'aliases', 'refseq']

    for index, row in issues_df.iterrows():

        gene = row['gene']
        acc = row['accession']
        acc_prefix = re.search('(.*?)(?=\.)', acc).group()

        try:
            hgnc_idx = hgnc_df.index[hgnc_df['refseq'] == acc_prefix]
            hgnc_row = hgnc_df.loc[hgnc_idx[0]]

            assert (gene in hgnc_row['approved']) or \
                (gene in hgnc_row['previous']) or \
                (gene in hgnc_row['aliases']), \
                f"{acc} wrongly matched to {hgnc_row['hgnc']}"

            hgnc = hgnc_row['hgnc']

        except IndexError:
            hgnc = '-'

        hgncs.append(hgnc)

    issues_df['hgnc'] = hgncs

    return issues_df


def unique_gene_list(issues_df):
    """ Generate a list of the unique genes with low coverage.

    args:
        issues_df [df]: data for exons with low coverage

    returns:
        problem_genes [list]
    """

    problem_genes = []

    for index, row in issues_df.iterrows():

        identifier = f"{row['gene']}\t{row['accession']}\t{row['hgnc']}"

        if identifier not in problem_genes:
            problem_genes.append(identifier)

    return problem_genes


def create_outputs(output_prefix, issues_df, problem_genes):
    """ Creates two outputs: a text file containing a list of unique
    genes with coverage issues, and an excel workbook detailing which
    exons of those genes were affected.

    args:
        output_file [str]: name for output xlsx workbook
        issues_df [df]: data for genes with coverage issues
        problem_genes [list]: unique genes with coverage issues
    """

    text_file = f"{output_prefix}.low_coverage_genes.txt"
    xlsx_file = f"{output_prefix}.low_coverage_gene_exons.xlsx"

    with open(text_file, 'w') as writer:
        writer.write(f"Sample: {output_prefix}\n\n")
        writer.write(f"Genes with <100% coverage at 30x:\n")
        writer.write("\n".join(problem_genes))

    issues_df.to_excel(xlsx_file, index=False)


def main():

    input_file = sys.argv[1]
    hgnc_dump = "230224_hgnc_dump.tsv"

    output_prefix, coverage_df = create_df(input_file)
    issues_df = find_coverage_issues(coverage_df)
    issues_df = get_hgncs(hgnc_dump, issues_df)
    problem_genes = unique_gene_list(issues_df)
    create_outputs(output_prefix, issues_df, problem_genes)


if __name__ == '__main__':
    main()
