from gnomad_hail.resources.resource_utils import TableResource
import hail as hl

lcr = TableResource(
    path="gs://gnomad-public/resources/intervals/LCRFromHengHg38.ht",
    import_func=hl.import_locus_intervals,
    import_args={
        "path": "gs://gnomad-public/resources/grch38/LCRFromHengHg38.txt",
        "reference_genome": 'GRCh38',
        "skip_invalid_intervals": True
    },
)

telomeres_and_centromeres = TableResource(
    path='gs://gnomad-public/resources/intervals/hg38.telomeresAndMergedCentromeres.ht',
    import_func=hl.import_bed,
    import_args={
        "path": 'gs://gnomad-public/resources/grch38/hg38.telomeresAndMergedCentromeres.bed',
        "reference_genome": 'GRCh38',
        "skip_invalid_intervals": True
    }
)
