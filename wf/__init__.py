"""latch/scrnaseq"""

import gzip
import subprocess
import sys
import types
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional, Tuple, Union

import lgenome
from dataclasses_json import dataclass_json
from flytekit import LaunchPlan
from latch import large_task, workflow, small_task
from latch.types import LatchDir, LatchFile

sys.stdout.reconfigure(line_buffering=True)


def run(cmd: List[str]):
    subprocess.run(cmd, check=True)


# TODO - patch latch with proper def __repr__ -> str
def ___repr__(self):
    return str(self.local_path)


LatchFile.__repr__ = types.MethodType(___repr__, LatchFile)


@dataclass_json
@dataclass
class SingleEndReads:
    r1: LatchFile


@dataclass_json
@dataclass
class PairedEndReads:
    r1: LatchFile
    r2: LatchFile


class ReadType(Enum):
    single = "single"
    paired = "paired"


class Strandedness(Enum):
    auto = "auto"


@dataclass_json
@dataclass
class Sample:
    name: str
    strandedness: Strandedness
    replicates: List[Union[SingleEndReads, PairedEndReads]]


class Technology(Enum):
    dropseq = "Drop-seq"
    chromiumv3 = "Chromium v3"
    chromiumv2 = "Chromium v2"
    gemcode = "GemCode"
    citeseq = "CiteSeq"
    celseq = "CelSeq"
    celseq2 = "CelSeq2"
    splitseqV1 = "SplitSeq v1"
    splitseqV2 = "SplitSeq v2"
    quartzseq2 = "QuartzSeq2"
    sciseq3 = "SciSeq3"


technology_to_flag = {
    Technology.dropseq: "--dropseq",
    Technology.chromiumv3: "--chromiumV3",
    Technology.chromiumv2: "--chromium",
    Technology.gemcode: "--gemcode",
    Technology.citeseq: "--citeseq",
    Technology.celseq: "--celseq",
    Technology.celseq2: "--celseq2",
    Technology.splitseqV1: "--splitseqV1",
    Technology.splitseqV2: "--splitseqV2",
    Technology.quartzseq2: "--quartzseq2",
    Technology.sciseq3: "--sciseq3",
}


class LatchGenome(Enum):
    RefSeq_hg38_p14 = "Homo sapiens (RefSeq hg38.p14)"
    RefSeq_T2T_CHM13v2_0 = "Homo sapiens (RefSeq T2T-CHM13v2.0)"
    RefSeq_R64 = "Saccharomyces cerevisiae (RefSeq R64)"
    RefSeq_GRCm39 = "Mus musculus (RefSeq GRCm39)"


# TODO - not used
@dataclass_json
@dataclass
class CustomGenome:
    gtf: LatchFile
    ref_genome: LatchFile
    ref_transcript: Optional[LatchFile]
    salmon_index: Optional[LatchFile]
    STAR_index: Optional[LatchFile]


class InsufficientCustomGenomeResources(Exception):
    pass


class MalformedSpliciIndex(Exception):
    pass


@small_task
def get_output_location(
    custom_output_dir: Optional[LatchDir],
    run_name: str
) -> str:
    run_name = run_name.rstrip("/").lstrip("/")
    if custom_output_dir is None:
        output_base = f"latch:///SCRNA-Seq Outputs/{run_name}/"
    else:
        remote_path: str = custom_output_dir.remote_path
        if remote_path[-1] != "/":
            remote_path += "/"
        output_base = f"{remote_path}{run_name}/"

    return output_base

@large_task
def make_splici_index(
    samples: List[Sample],
    read_length: Optional[int],
    splici_index: Optional[LatchDir],
    latch_genome: LatchGenome,
    output_name: str,
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
) -> Tuple[LatchDir, LatchFile]:

    if splici_index is not None:
        print("Using provided splici index")
        for file in Path(splici_index).iterdir():
            if file.name.endswith("_t2g_3col.tsv"):
                return splici_index, file
        raise MalformedSpliciIndex(
            "Could not find ...t2g_3col.tsv in provided splici index"
        )

    if read_length is None:
        print("No read length provided.\n\tEstimating read length...")
        reads_file = Path(samples[0].replicates[0].r2)
        lens = []
        if reads_file.suffix == ".gz":
            with gzip.open(reads_file, "r") as f:
                for i, l in enumerate(f):
                    if i > 400:
                        break
                    if i % 4 != 1:
                        continue
                    lens.append(len(l))
        else:
            with open(reads_file, "r") as f:
                for i, l in enumerate(f):
                    if i > 400:
                        break
                    if i % 4 != 1:
                        continue
                    lens.append(len(l))

        read_length = int(sum(lens) / len(lens))

    print(f"\tRead Length: {read_length}")

    if custom_gtf is not None and custom_ref_genome is not None:
        print("Using custom genome and GTF file to generate splici index")

        genome_path = Path(custom_ref_genome)
        gtf_path = Path(custom_gtf)
    else:
        print("Using Latch genome to generate splici index")
        gm = lgenome.GenomeManager(latch_genome.name)
        genome_path = gm.download_ref_genome()
        gtf_path = gm.download_gtf()

    print("Generating splici transcriptome...")

    flank_trim_length = 5

    print(f"\tUsing flank trim length: {flank_trim_length}")

    txome_cmd = [
        "pyroe",
        "make-splici",
        str(genome_path),
        str(gtf_path),
        str(read_length),
        "splici_txome",
        "--flank-trim-length",
        f"{flank_trim_length}",
        "--filename-prefix",
        "splici_txome",
        "--dedup-seqs",
    ]

    txome_process = subprocess.Popen(
        " ".join(txome_cmd),
        shell=True,
        errors="replace",
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    while True:
        realtime_output = txome_process.stdout.readline()

        if realtime_output == "" and txome_process.poll() is not None:
            break
        if realtime_output:
            print("\t" + realtime_output.strip())

    retval = txome_process.poll()
    if retval != 0:
        raise Exception(f"\tPyroe failed with exit code {retval}")

    print("Splici transcriptome generated.")

    print("Generating splici index...")

    index_cmd = [
        "salmon",
        "index",
        "-t",
        f"splici_txome/splici_txome_fl{read_length - flank_trim_length}.fa",
        "-i",
        "splici_index",
        "-p",
        "96",
    ]

    index_process = subprocess.Popen(
        " ".join(index_cmd),
        shell=True,
        errors="replace",
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    while True:
        realtime_output = index_process.stdout.readline()

        if realtime_output == "" and index_process.poll() is not None:
            break
        if realtime_output:
            print("\t" + realtime_output.strip())

    retval = index_process.poll()
    if retval != 0:
        raise Exception(f"\tSalmon index failed with exit code {retval}")

    print("Splici index generated. Packaging Files...")
    return LatchDir(
        "/root/splici_index", f"{output_name}splici_index"
    ), LatchFile(
        f"splici_txome/splici_txome_fl{read_length - flank_trim_length}_t2g_3col.tsv",
        f"{output_name}splici_index/splici_txome_fl{read_length - flank_trim_length}_t2g_3col.tsv",
    )


@large_task
def map_reads(
    splici_index: LatchDir,
    samples: List[Sample],
    technology: Technology,
    output_name: str,
) -> LatchDir:
    sample_args = []
    if type(samples[0].replicates[0]) is SingleEndReads:
        sample_args.append("-r")
        for sample in samples:
            for replicate in sample.replicates:
                sample_args.append(f"{Path(replicate.r1)}")
    else:
        r1 = []
        r2 = []
        for sample in samples:
            for replicate in sample.replicates:
                r1.append(str(Path(replicate.r1)))
                r2.append(str(Path(replicate.r2)))
        sample_args.append("-1")
        sample_args += r1
        sample_args.append("-2")
        sample_args += r2

    # we allow the library type to be inferred via `-l A` flag.
    alevin_cmd = [
        "salmon",
        "alevin",
        "-p",
        "96",
        "-i",
        str(Path(splici_index)),
        "-l",
        "A",
        technology_to_flag.get(technology),
        "--rad",
    ]
    alevin_cmd += sample_args

    alevin_cmd += [
        "-o",
        "map",
    ]

    print("Mapping reads to transcripts...")

    alevin_process = subprocess.Popen(
        " ".join(alevin_cmd),
        shell=True,
        errors="replace",
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    while True:
        realtime_output = alevin_process.stdout.readline()

        if realtime_output == "" and alevin_process.poll() is not None:
            break
        if realtime_output:
            print("\t" + realtime_output.strip())

    retval = alevin_process.poll()
    if retval != 0:
        raise Exception(f"\tSalmon alevin failed with exit code {retval}")

    print("Transcript mapping complete. Packaging Files...")
    return LatchDir("/root/map", f"{output_name}intermediate_mapping")


@large_task
def quantify_reads(
    map_dir: LatchDir,
    tg_map: LatchFile,
    output_name: str,
) -> LatchDir:
    print("Generating permit list...")
    permit_list_cmd = [
        "alevin-fry",
        "generate-permit-list",
        "-d",
        "fw",
        "-k",
        "-i",
        str(Path(map_dir)),
        "-o",
        "quant",
    ]

    permit_list_process = subprocess.Popen(
        " ".join(permit_list_cmd),
        shell=True,
        errors="replace",
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    while True:
        realtime_output = permit_list_process.stdout.readline()

        if realtime_output == "" and permit_list_process.poll() is not None:
            break
        if realtime_output:
            print("\t" + realtime_output.strip())

    retval = permit_list_process.poll()
    if retval != 0:
        raise Exception(
            f"\tAlevin-fry generate-permit-list failed with exit code {retval}"
        )

    print("Collating RAD files...")
    collate_cmd = [
        "alevin-fry",
        "collate",
        "-t",
        "96",
        "-i",
        "quant",
        "-r",
        str(Path(map_dir)),
    ]

    collate_process = subprocess.Popen(
        " ".join(collate_cmd),
        shell=True,
        errors="replace",
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    while True:
        realtime_output = collate_process.stdout.readline()

        if realtime_output == "" and collate_process.poll() is not None:
            break
        if realtime_output:
            print("\t" + realtime_output.strip())

    retval = collate_process.poll()
    if retval != 0:
        raise Exception(f"\tAlevin-fry collate failed with exit code {retval}")

    print("Quantifying reads...")
    quant_cmd = [
        "alevin-fry",
        "quant",
        "-t",
        "96",
        "-i",
        "quant",
        "-o",
        "quant_res",
        "--tg-map",
        str(Path(tg_map)),
        "--resolution",
        "cr-like",
        "--use-mtx",
    ]

    quant_process = subprocess.Popen(
        " ".join(quant_cmd),
        shell=True,
        errors="replace",
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    while True:
        realtime_output = quant_process.stdout.readline()

        if realtime_output == "" and quant_process.poll() is not None:
            break
        if realtime_output:
            print("\t" + realtime_output.strip())

    retval = quant_process.poll()
    if retval != 0:
        raise Exception(f"\tAlevin-fry quant failed with exit code {retval}")

    print("Quantification complete. Packaging Files...")
    return LatchDir("/root/quant_res", f"{output_name}counts")


@workflow
def scrnaseq(
    samples: List[Sample],
    read_length: Optional[int],
    run_name: str,
    ta_ref_genome_fork: str,
    latch_genome: LatchGenome,
    technology: Technology = Technology.chromiumv3,
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
    splici_index: Optional[LatchDir] = None,
    custom_output_dir: Optional[LatchDir] = None,
) -> LatchDir:
    """Performs alignment & quantification on Single Cell RNA-Sequencing reads.

    Single Cell RNA-Seq (Alignment & Quantification)
    ----

    This workflow will produce gene and transcript counts from single cell RNA-seq
    sample reads.

    This current iteration of the workflow has three steps:

    1. Generating an intron aware index [Splici Index](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2.supplementary-material)
    2. Selective alignment and quantification of reads using [Salmon](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2)
    3. Quality Control using [alevinQC](https://github.com/csoneson/alevinQC)

    ### More about Selective Alignment

    This workflow uses Salmon's selective alignment described in this
    [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8),
    which achieves greater accuracy than traditional alignment methods while
    using fewer computational resources.


    __metadata__:
        display_name: Single Cell RNAseq
        author:
            name: LatchBio
            email:
            github:
        repository: github.com/latchbio/scrnaseq
        license:
            id: MIT
        flow:
        - section: Samples
          flow:
            - text: >-
                  Sample files can be provided and their read type can be
                  inferred from their name or this information can be specified manually.
                  Sample strandedness is inferred automatically (learn more).

            - params:
                - samples
                - read_length
                - technology
        - section: Alignment & Quantification
          flow:
            - text: >-
                  This workflow uses Salmon's selective alignment described in this
                  [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8),
                  which achieves greater accuracy than traditional alignment methods while
                  using less computational resources. To fix intron problem...

            - fork: alignment_quantification_tools
              flows:
                traditional:
                    display_name: Selective Alignment
                    flow:
                        - fork: ta_ref_genome_fork
                          flows:
                            database:
                                display_name: Curated Genome For Index Generation
                                flow:
                                    - text: >-
                                        We have curated a set of reference
                                        genome data for ease and
                                        reproducibility. More information about
                                        these managed files can be found
                                        [here](https://github.com/latchbio/latch-genomes).
                                    - params:
                                        - latch_genome
                            custom:
                                display_name: Custom Genome For Index Generation
                                _tmp_unwrap_optionals:
                                    - custom_gtf
                                    - custom_ref_genome
                                flow:
                                    - params:
                                        - custom_ref_genome
                                        - custom_gtf
                            prebuilt:
                                display_name: Prebuilt Splici Index
                                _tmp_unwrap_optionals:
                                    - splici_index
                                flow:
                                    - params:
                                        - splici_index

        - section: Output Settings
          flow:
          - params:
              - run_name
          - fork: output_location_fork
            flows:
                default:
                    display_name: Default
                    flow:
                    - text:
                        Output will be at default location in the data
                        viewer - RNA-Seq Outputs/"Run Name"
                custom:
                    display_name: Specify Custom Path
                    _tmp_unwrap_optionals:
                        - custom_output_dir
                    flow:
                    - params:
                        - custom_output_dir
    Args:

        samples:
            Here you can organize your FastQ files by sample and add technical
            replicates for each sample.  Biological replicates should be
            organized as separate samples.

          __metadata__:
            display_name: Sample Sheet
            batch_table_column: true
            _tmp:
                custom_ingestion: auto

        read_length:
            The length of the reads in the sample. This is used to generate the splici index.
            If you are providing an index, you can ignore this parameter. If you do not provide
            an index or fill out this parameter, the read length will be estimated from the
            first 100 reads in the first sample.

            __metadata__:
                display_name: Reads Length
                batch_table_column: true


        alignment_quantification_tools:

          __metadata__:
            display_name: Alignment & Quantification Method

        ta_ref_genome_fork:

            __metadata__:
                display_name: Reference Genome Source

        latch_genome:
          Curated reference files for specific genome sources and builds.

          __metadata__:
            batch_table_column: true
            display_name: Genome Database Option

        technology:
            Sequencing technology used for this experiment. If the technology you used
            is not an option, reach out of aidan@latch.bio so we can implement it.

            __metadata__:
                display_name: Sequencing Technology

        custom_ref_genome:
          The reference genome you want to align you samples to.

          __metadata__:
            display_name: Reference Genome File
            appearance:
                detail: (.fasta, .fasta.gz, .fa, .fa.gz, .fna, .fna.gz)

        custom_gtf:
          The gene annonation file that corresponds to the reference genome
          provided.

          __metadata__:
            display_name: Annotation File
            appearance:
                detail: (.gtf)

        splici_index:
            You are able to provide a prebuilt splici index directory for your genome.
            This will speed up run time as the index is generated if none is provided.
            This should be the output from a previous scrnaseq run or generated using
            this [tutorial](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/).

          __metadata__:
            display_name: Splici index

        save_indices:
            If you provided a custom genome you can output the alignment
            indexes generated from this run for use in future runs. This will
            speed up runtime since the workflow doesn't have to then regenerate
            the indexes.

          __metadata__:
            display_name: Save Generated Reference Indexes

        run_name:
          A name for this analysis run, this will be used to name outputs from
          this run.

          __metadata__:
            batch_table_column: true
            display_name: Run Name

        output_location_fork:

        custom_output_dir:
          You can provide a custom location where this run's analysis outputs
          will be located.

          __metadata__:
            display_name: Custom Output Location
    """

    output_name = get_output_location(
        custom_output_dir=custom_output_dir,
        run_name=run_name,
    )

    (splici_index, t2g) = make_splici_index(
        samples=samples,
        read_length=read_length,
        splici_index=splici_index,
        latch_genome=latch_genome,
        custom_gtf=custom_gtf,
        custom_ref_genome=custom_ref_genome,
        output_name=output_name,
    )

    mapped_reads = map_reads(
        splici_index=splici_index,
        samples=samples,
        technology=technology,
        output_name=output_name,
    )

    quantified_reads = quantify_reads(
        map_dir=mapped_reads,
        tg_map=t2g,
        output_name=output_name,
    )

    return quantified_reads


if __name__ == "wf":
    LaunchPlan.create(
        "wf.__init__.scrnaseq.Fast Test Data",
        scrnaseq,
        default_inputs={
            "samples": [
                Sample(
                    name="dummy_data",
                    replicates=[
                        PairedEndReads(
                            r1=LatchFile(
                                "s3://latch-public/welcome/scrnaseq/dummy_L001_R1_001.fastq",
                            ),
                            r2=LatchFile(
                                "s3://latch-public/welcome/scrnaseq/dummy_L001_R2_001.fastq",
                            ),
                        ),
                        PairedEndReads(
                            r1=LatchFile(
                                "s3://latch-public/welcome/scrnaseq/dummy_L002_R1_001.fastq",
                            ),
                            r2=LatchFile(
                                "s3://latch-public/welcome/scrnaseq/dummy_L002_R2_001.fastq",
                            ),
                        ),
                    ],
                    strandedness=Strandedness.auto,
                )
            ],
            "run_name": "Fast Test Data",
            "custom_gtf": LatchFile(
                "s3://latch-public/welcome/scrnaseq/dummy_genes.gtf"
            ),
            "custom_ref_genome": LatchFile(
                "s3://latch-public/welcome/scrnaseq/dummy_genome.fa"
            ),
        },
    )
