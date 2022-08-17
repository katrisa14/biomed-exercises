from typing import List
import click
from .normalize_vcf import normalize_vcf


@click.group()
def cli():
    pass


@cli.command()
@click.argument("vcf_in", required=True)
@click.argument("vcf_out", required=True)
def normalize(vcf_in: str, vcf_out: str):
    normalize_vcf(vcf_in, vcf_out)
