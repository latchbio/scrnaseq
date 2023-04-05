"""latch/scrnaseq"""

import gzip
import shutil
import subprocess
import sys
import time
import traceback
import types
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from shlex import quote
from typing import Callable, Dict, List, Optional, TextIO, Tuple, Union

import anndata
import deepsort
import lgenome
import mygene
import pyroe
import scanpy as sc
from dataclasses_json import dataclass_json
from flytekit import LaunchPlan
from latch import large_task, message, small_task, workflow
from latch.types import LatchDir, LatchFile

sys.stdout.reconfigure(line_buffering=True)


def is_gzipped(path: Path) -> bool:
    return path.suffix == ".gz"


def safe_open_fxn(path: Path) -> Callable[[Path, str], TextIO]:
    return gzip.open if is_gzipped(path) else open


# TODO - patch latch with proper def __repr__ -> str
def ___repr__(self):
    return str(self.local_path)


LatchFile.__repr__ = types.MethodType(___repr__, LatchFile)


MYGENE_SCOPES = [
    "entrezgene",
    "ensembl.gene",
    "symbol",
    "name",
    "alias",
    "refseq",
    "unigene",
    "homologene",
    "accession",
    "ensembl.transcript",
    "ensembl.protein",
    "uniprot",
    "pdb",
    "prosite",
    "pfam",
    "interpro",
    "mim",
    "pharmgkb",
    "reporter",
    "reagent",
    "go",
    "hgnc",
    "hprd",
    "mgi",
    "rgd",
    "flybase",
    "wormbase",
    "zfin",
    "tair",
    "xenbase",
    "mirbase",
    "retired",
]


class Species(Enum):
    human = "Human"
    mouse = "Mouse"


class HumanTissue(Enum):
    adipose = "adipose"
    adrenal_gland = "adrenal gland"
    artery = "artery"
    ascending_colon = "ascending colon"
    bladder = "bladder"
    blood = "blood"
    bone_marrow = "bone marrow"
    brain = "brain"
    cervix = "cervix"
    chorionic_villus = "chorionic villus"
    colorectum = "colorectum"
    cord_blood = "cord blood"
    epityphlon = "epityphlon"
    esophagus = "esophagus"
    fallopian_tube = "fallopian tube"
    female_gonad = "female gonad"
    fetal_adrenal_gland = "fetal adrenal gland"
    fetal_brain = "fetal brain"
    fetal_calvaria = "fetal calvaria"
    fetal_eye = "fetal eye"
    fetal_heart = "fetal heart"
    fetal_intestine = "fetal intestine"
    fetal_kidney = "fetal kidney"
    fetal_liver = "fetal liver"
    fetal_lung = "fetal lung"
    fetal_male_gonad = "fetal male gonad"
    fetal_muscle = "fetal muscle"
    fetal_pancreas = "fetal pancreas"
    fetal_rib = "fetal rib"
    fetal_skin = "fetal skin"
    fetal_spinal_cord = "fetal spinal cord"
    fetal_stomach = "fetal stomach"
    fetal_thymus = "fetal thymus"
    gall_bladder = "gall bladder"
    heart = "heart"
    kidney = "kidney"
    liver = "liver"
    lung = "lung"
    muscle = "muscle"
    neonatal_adrenal_gland = "neonatal adrenal gland"
    omentum = "omentum"
    pancreas = "pancreas"
    placenta = "placenta"
    pleura = "pleura"
    prostat = "prostat"
    spleen = "spleen"
    stomach = "stomach"
    temporal_lobe = "temporal lobe"
    thyroid = "thyroid"
    trachea = "trachea"
    ureter = "ureter"


class MouseTissue(Enum):
    bladder = "bladder"
    blood = "blood"
    bone_marrow = "bone marrow"
    bone_marrow_mesenchyme = "bone marrow mesenchyme"
    brain = "brain"
    embryonic_mesenchyme = "embryonic mesenchyme"
    fetal_brain = "fetal brain"
    fetal_intestine = "fetal intestine"
    fetal_liver = "fetal liver"
    fetal_lung = "fetal lung"
    fetal_stomach = "fetal stomach"
    intestine = "intestine"
    kidney = "kidney"
    liver = "liver"
    lung = "lung"
    mammary_gland = "mammary gland"
    muscle = "muscle"
    neonatal_calvaria = "neonatal calvaria"
    neonatal_heart = "neonatal heart"
    neonatal_muscle = "neonatal muscle"
    neonatal_pancreas = "neonatal pancreas"
    neonatal_rib = "neonatal rib"
    neonatal_skin = "neonatal skin"
    ovary = "ovary"
    pancreas = "pancreas"
    placenta = "placenta"
    prostate = "prostate"
    spleen = "spleen"
    stomach = "stomach"
    testis = "testis"
    thymus = "thymus"
    uterus = "uterus"


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


@dataclass
class ProtocolGeometry:
    barcode_length: int
    umi_length: int


technology_to_geometry = {
    Technology.dropseq: ProtocolGeometry(barcode_length=12, umi_length=8),
    Technology.chromiumv3: ProtocolGeometry(barcode_length=16, umi_length=12),
    Technology.chromiumv2: ProtocolGeometry(barcode_length=16, umi_length=10),
    Technology.gemcode: ProtocolGeometry(barcode_length=14, umi_length=10),
    Technology.citeseq: ProtocolGeometry(barcode_length=16, umi_length=10),
    Technology.celseq: ProtocolGeometry(barcode_length=8, umi_length=6),
    Technology.celseq2: ProtocolGeometry(barcode_length=6, umi_length=6),
    Technology.splitseqV1: ProtocolGeometry(barcode_length=24, umi_length=10),
    Technology.splitseqV2: ProtocolGeometry(barcode_length=24, umi_length=10),
    Technology.quartzseq2: ProtocolGeometry(barcode_length=15, umi_length=8),
    Technology.sciseq3: ProtocolGeometry(barcode_length=21, umi_length=8),
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


class TranscriptomeGenerationError(Exception):
    pass


class GTFGenerationError(Exception):
    pass


class IndexGenerationError(Exception):
    pass


class SalmonAlevinError(Exception):
    pass


class AlevinFryQuantError(Exception):
    pass


class AlevinFryCollateError(Exception):
    pass


class AlevinFryGeneratePermitListError(Exception):
    pass


class H5ADGenerationError(Exception):
    pass


class ReportGenerationError(Exception):
    pass


class SingleEndReadsUnsupported(Exception):
    pass


class SampleToCellularBarcodeMapError(Exception):
    pass


@small_task
def get_output_location(custom_output_dir: Optional[LatchDir], run_name: str) -> str:
    print("Generating output location...")
    stripped_run_name = run_name.rstrip("/").lstrip("/")

    if stripped_run_name == "":
        raise RuntimeError(f"Run name must be provided and contain alphanumerical characters: found `{run_name}`")

    if custom_output_dir is None:
        output_base = f"latch:///SCRNA-Seq Outputs/{stripped_run_name}/"
    else:
        remote_path: str = custom_output_dir.remote_path
        if remote_path[-1] != "/":
            remote_path += "/"
        output_base = f"{remote_path}{stripped_run_name}/"

    print("\tOutput location: " + output_base)
    message("info", {"title": "Output Directory", "body": output_base})
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
        message(
            "info",
            {
                "title": "Using provided splici index",
                "body": "Skipping index generation",
            },
        )
        for file in Path(splici_index).iterdir():
            if file.name.endswith("_t2g_3col.tsv"):
                return splici_index, file
        message(
            "error",
            {
                "title": "Malformed Splici Index",
                "body": "Could not find ...t2g_3col.tsv in provided splici index",
            },
        )
        raise MalformedSpliciIndex(
            "Could not find ...t2g_3col.tsv in provided splici index"
        )

    print("Estimating read length...")
    if read_length is None:
        if isinstance(samples[0].replicates[0], SingleEndReads):
            reads_file = Path(samples[0].replicates[0].r1)
        else:
            reads_file = Path(samples[0].replicates[0].r2)
        lens = []
        with safe_open_fxn(reads_file)(reads_file, "rt") as f:
            for i, l in enumerate(f):
                if i > 400:
                    break
                if i % 4 != 1:
                    continue
                lens.append(len(l))

        read_length = int(sum(lens) / len(lens))

    message("info", {"title": "Computed Read Length", "body": f"{read_length}"})
    print(f"\tRead Length: {read_length}")

    if custom_gtf is not None and custom_ref_genome is not None:
        message(
            "info",
            {
                "title": "Generating Index Using Custom Genome",
                "body": f"Genome: {custom_ref_genome.remote_path}\nGTF: {custom_gtf.remote_path}",
            },
        )
        print("Using custom genome and GTF file to generate splici index")

        genome_path = Path(custom_ref_genome)
        gtf_path = Path(custom_gtf)

        if is_gzipped(genome_path):
            print("\tUnzipping genome file")
            message("info", {"title": "Unzipping Genome", "body": ""})
            with gzip.open(genome_path, "rb") as f_in:
                with open(genome_path.with_suffix(""), "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            genome_path = genome_path.with_suffix("")

        if is_gzipped(gtf_path):
            print("\tUnzipping GTF file")
            message("info", {"title": "Unzipping GTF", "body": ""})
            with gzip.open(gtf_path, "rb") as f_in:
                with open(gtf_path.with_suffix(""), "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            gtf_path = gtf_path.with_suffix("")

        if (
            gtf_path.suffix == ".gff3"
            or gtf_path.suffix == ".gff"
            or gtf_path.suffix == ".gff2"
        ):
            print("\tConverting GFF to GTF3")
            message(
                "warning",
                {"title": "Attempting to convert GFF to GTF", "body": f"{gtf_path}"},
            )
            new_gtf_path = gtf_path.with_suffix(".gtf")

            gffread_cmd = [
                "gffread",
                str(gtf_path),
                "-T",
                "-o",
                str(new_gtf_path),
            ]

            gffread_cmd = " ".join(gffread_cmd)
            print("Command: " + gffread_cmd)
            message("info", {"title": "Command", "body": gffread_cmd})

            gffread_process = subprocess.Popen(
                gffread_cmd,
                shell=True,
                errors="replace",
                encoding="utf-8",
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            while True:
                realtime_output = gffread_process.stdout.readline()

                if realtime_output == "" and gffread_process.poll() is not None:
                    break
                if realtime_output:
                    print("\t" + realtime_output.strip())

            retval = gffread_process.poll()
            if retval != 0:
                message(
                    "error",
                    {
                        "title": "GTF Generation Error",
                        "body": f"View logs to see error",
                    },
                )
                raise GTFGenerationError(f"\tGffread failed with error code {retval}")

            gtf_path = new_gtf_path

    else:
        message(
            "info",
            {
                "title": "Generating Index Using Latch Genome",
                "body": f"Name: {latch_genome.name}",
            },
        )
        print("Using Latch genome to generate splici index")
        gm = lgenome.GenomeManager(latch_genome.name)
        genome_path = gm.download_ref_genome()
        gtf_path = gm.download_gtf()

    print("Generating splici transcriptome...")

    flank_trim_length = 5

    print(f"\tUsing flank trim length: {flank_trim_length}")
    message("info", {"title": "Flank Trim Length", "body": f"{flank_trim_length}"})

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

    txome_cmd = " ".join(txome_cmd)
    print("Command: " + txome_cmd)
    message("info", {"title": "Command", "body": txome_cmd})

    pyroe_attempts = 0
    while True:
        try:
            # linear backoff of failed network requests
            time.sleep(pyroe_attempts)
            print("Pyroe Attempt: " + str(pyroe_attempts + 1))
            txome_process = subprocess.Popen(
                txome_cmd,
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
            if retval == 0:
                break
            else:
                raise TranscriptomeGenerationError(
                    f"\tPyroe failed with exit code {retval}"
                )
        except TranscriptomeGenerationError as e:
            pyroe_attempts += 1
            if pyroe_attempts == 5:
                message(
                    "error",
                    {
                        "title": "Transcriptome Generation Error",
                        "body": f"View logs to see error",
                    },
                )
                raise e

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

    index_cmd = " ".join(index_cmd)
    print("Command: " + index_cmd)
    message("info", {"title": "Command", "body": index_cmd})

    index_process = subprocess.Popen(
        index_cmd,
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
        message(
            "error",
            {"title": "Index Generation Error", "body": f"View logs to see error"},
        )
        raise IndexGenerationError(f"\tSalmon index failed with exit code {retval}")

    message("info", {"title": "Success", "body": ""})
    print("Splici index generated. Packaging Files...")
    return (
        LatchDir("/root/splici_index", f"{output_name}splici_index"),
        LatchFile(
            f"splici_txome/splici_txome_fl{read_length - flank_trim_length}_t2g_3col.tsv",
            f"{output_name}splici_index/splici_txome_fl{read_length - flank_trim_length}_t2g_3col.tsv",
        ),
    )


@large_task
def map_reads(
    splici_index: LatchDir,
    samples: List[Sample],
    technology: Technology,
    output_name: str,
) -> Tuple[List[LatchDir], List[str]]:
    nrof_samples = 0

    mapping_dirs: List[LatchDir] = []
    sample_names: List[str] = []
    for sample in samples:
        nrof_samples += 1
        sample_names.append(sample.name)
        message(
            "info",
            {
                "title": f"Downloading Sample Data For {sample.name}",
                "body": "",
            },
        )
        print(f"Downloading sample data...")
        if isinstance(sample.replicates[0], SingleEndReads):
            raise SingleEndReadsUnsupported("Single end reads are not supported")

        r1 = []
        r2 = []
        sample_args = []
        for replicate in sample.replicates:
            r1.append(str(Path(replicate.r1)))
            r2.append(str(Path(replicate.r2)))
        sample_args.append("-1")
        sample_args += r1
        sample_args.append("-2")
        sample_args += r2

        message(
            "info",
            {
                "title": f"Success",
                "body": "",
            },
        )
        print("Done")

        message(
            "info",
            {
                "title": f"Mapping {sample.name} To Transcripts",
                "body": "",
            },
        )
        print(f"Mapping {sample.name} to transcripts...")

        # we allow the library type to be inferred via `-l A` flag.
        alevin_cmd = [
            "salmon",
            "alevin",
            "-p",
            "96",
            "-i",
            str(Path(splici_index)),
            "-lA",
            technology_to_flag.get(technology),
            "--rad",
        ]
        alevin_cmd += sample_args

        alevin_cmd += [
            "-o",
            f"{sample.name}_map",
        ]

        alevin_cmd = " ".join(alevin_cmd)
        print("Command: " + alevin_cmd)
        message("info", {"title": "Command", "body": alevin_cmd})

        alevin_process = subprocess.Popen(
            alevin_cmd,
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
            message(
                "error",
                {"title": "Salmon Alevin Error", "body": f"View logs to see error"},
            )
            raise SalmonAlevinError(f"\tSalmon alevin failed with exit code {retval}")

        mapping_dirs.append(
            LatchDir(
                f"/root/{sample.name}_map",
                f"{output_name}{sample.name}/intermediate_mapping",
            )
        )
        message("info", {"title": "Success", "body": ""})

    message("info", {"title": "Success", "body": f"{nrof_samples} samples mapped"})
    print(f"Transcript mapping complete for {nrof_samples}. Packaging Files...")
    return mapping_dirs, sample_names


@large_task
def quantify_reads(
    mapping_dirs: List[LatchDir],
    sample_names: List[str],
    tg_map: LatchFile,
    output_name: str,
) -> Tuple[List[LatchDir], List[LatchDir]]:
    preprocessing_dirs: List[LatchDir] = []
    quantification_dirs: List[LatchDir] = []
    for sample_name, mapping_dir in zip(sample_names, mapping_dirs):
        message("info", {"title": f"Quantifying {sample_name}", "body": ""})
        message(
            "info",
            {"title": f"Generating Permit List For Sample {sample_name}", "body": ""},
        )
        print(f"Generating permit list...")
        permit_list_cmd = [
            "alevin-fry",
            "generate-permit-list",
            "-d",
            "either",
            "-k",
            "-i",
            str(Path(mapping_dir)),
            "-o",
            f"{sample_name}_quant",
        ]

        permit_list_cmd = " ".join(permit_list_cmd)
        print("Command: " + permit_list_cmd)
        message("info", {"title": "Command", "body": permit_list_cmd})

        permit_list_process = subprocess.Popen(
            permit_list_cmd,
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
            message(
                "error",
                {
                    "title": "AlevinFry Generate Permit List Error",
                    "body": f"View logs to see error",
                },
            )
            raise AlevinFryGeneratePermitListError(
                f"\tAlevin-fry generate-permit-list failed with exit code {retval} for {sample_name}"
            )
        message("info", {"title": "Success", "body": ""})

        message(
            "info",
            {"title": f"Collating RAD File For Sample {sample_name}", "body": ""},
        )
        print("Collating RAD files...")
        collate_cmd = [
            "alevin-fry",
            "collate",
            "-t",
            "96",
            "-i",
            f"{sample_name}_quant",
            "-r",
            str(Path(mapping_dir)),
        ]

        collate_cmd = " ".join(collate_cmd)
        print("Command: " + collate_cmd)
        message("info", {"title": "Command", "body": collate_cmd})

        collate_process = subprocess.Popen(
            collate_cmd,
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
            message(
                "error",
                {"title": "AlevinFry Collate Error", "body": f"View logs to see error"},
            )
            raise AlevinFryCollateError(
                f"\tAlevin-fry collate failed with exit code {retval} for {sample_name}"
            )

        message("info", {"title": "Success", "body": ""})

        message(
            "info", {"title": f"Quantifying Reads For Sample {sample_name}", "body": ""}
        )
        print("Quantifying reads...")
        quant_cmd = [
            "alevin-fry",
            "quant",
            "--num-bootstraps",
            "0",
            "-t",
            "96",
            "-i",
            f"{sample_name}_quant",
            "-o",
            f"{sample_name}_quant_res",
            "--tg-map",
            str(Path(tg_map)),
            "--resolution",
            "cr-like",
            "--use-mtx",
            "--summary-stat",
        ]

        quant_cmd = " ".join(quant_cmd)
        print("Command: " + quant_cmd)
        message("info", {"title": "Command", "body": quant_cmd})

        quant_process = subprocess.Popen(
            quant_cmd,
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
            message(
                "error", {"title": "AlevinFry Quant", "body": f"View logs to see error"}
            )
            raise AlevinFryQuantError(
                f"\tAlevin-fry quant failed with exit code {retval}"
            )

        preprocessing_dirs.append(
            LatchDir(
                f"/root/{sample_name}_quant",
                f"{output_name}{sample_name}/quant_preprocessing",
            )
        )
        quantification_dirs.append(
            LatchDir(
                f"/root/{sample_name}_quant_res",
                f"{output_name}{sample_name}/raw_counts",
            )
        )

        message("info", {"title": "Success", "body": f"{sample_name}"})

    print("Quantification complete. Packaging Files...")
    return preprocessing_dirs, quantification_dirs


def infer_cell_type(
    h5ad: anndata.AnnData,
    mouse_tissue: Optional[MouseTissue],
    human_tissue: Optional[HumanTissue],
) -> Optional[List[Tuple[str, str]]]:
    """
    Infer the cell type if mouse tissue or human tissue selection.

    Returns a list of tuples of the form (cell_type, cell_subtype)
    """
    # Infer the cell type
    species = None
    if mouse_tissue:
        species = "mouse"
        tissue = "_".join(mouse_tissue.value.split(" ")).capitalize()
    elif human_tissue:
        species = "human"
        tissue = "_".join(human_tissue.value.split(" ")).capitalize()
    else:
        return None

    message("info", {"title": f"Running Cell Type Analysis", "body": ""})
    try:
        sc.pp.normalize_per_cell(h5ad)
        h5ad = h5ad.T
        h5ad.write_csvs("cell_type.csv", skip_data=False)

        cell_names = []
        with open("cell_type/var.csv") as f:
            f.readline()
            for line in f.readlines():
                cell_names.append(line.split(",")[0])

        with open("cell_type/inferme.csv", "w") as out:
            out.write("," + ",".join(cell_names) + "\n")
            with open("cell_type/X.csv") as counts:
                with open("cell_type/obs.csv") as symbols:
                    symbols.readline()
                    for symbol, raw in zip(symbols.readlines(), counts.readlines()):
                        symbol = symbol.split(",")[1].strip()
                        if symbol == "NA":
                            continue
                        out.write((",".join([symbol, raw])).strip() + "\n")

        model = deepsort.DeepSortPredictor(species=species, tissue=tissue)

        print("Running DeepSort...")
        model.predict("cell_type/inferme.csv", save_path="cell_type_results")

        with open(f"cell_type_results/{species}_{tissue}_inferme.csv") as f:
            cell_types = [
                (line.split(",")[1].strip(), line.split(",")[2].strip())
                for line in f.readlines()
            ]

        shutil.rmtree("/root/cell_type", ignore_errors=True)
        shutil.rmtree("/root/cell_type_results", ignore_errors=True)

        message("info", {"title": f"Done", "body": ""})
        return cell_types[1:]
    except Exception as e:
        message("warning", {"title": f"Cell Type Analysis Error", "body": str(e)})
        traceback.print_exc()
        return None


def run_decontX(
    h5ad_path: Path,
    sample_name: str,
) -> None:
    """
    Remove ambient RNA signal from count matrix using decontX.
    Executes in place on an input file.
    """
    message("info", {"title": f"Running decontX on {sample_name}", "body": ""})
    print(f"Running decontX on {sample_name}...")
    decontx_cmd = [
        "Rscript",
        "decontx.R",
        str(h5ad_path),
    ]

    decontx_process = subprocess.Popen(
        " ".join(decontx_cmd),
        shell=True,
        errors="replace",
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    while True:
        realtime_output = decontx_process.stdout.readline()

        if realtime_output == "" and decontx_process.poll() is not None:
            break
        if realtime_output:
            print("\t" + realtime_output.strip())

    retval = decontx_process.poll()
    if retval != 0:
        message(
            "warning",
            {"title": "DecontX Ambient RNA Removal", "body": f"View logs to see error"},
        )
        print(f"\tDecontX failed with exit code {retval}")


@large_task
def h5ad(
    quant_dirs: List[LatchDir],
    output_name: str,
    sample_names: List[str],
    mouse_tissue: Optional[MouseTissue],
    human_tissue: Optional[HumanTissue],
    remove_ambient_rna: bool = False,
) -> Tuple[List[LatchFile], LatchFile]:
    message("info", {"title": f"Mining For Gene Metadata", "body": ""})
    print("Mining for gene metadata...")

    h5ad_output_standard: anndata.AnnData = pyroe.load_fry(
        frydir=str(Path(quant_dirs[0])), output_format="scRNA"
    )
    gene_names = h5ad_output_standard.var_names

    mg = mygene.MyGeneInfo()
    res = mg.querymany(gene_names, scopes=",".join(MYGENE_SCOPES), fields="all")

    gene_info_by_input = {}
    for result in res:
        try:
            gene_info_by_input[result["query"]] = result
        except:
            print(f"Odd result in gene mining: {result}")
            pass

    gene_symbols = []
    gene_types = []
    gene_long_names = []
    found_symbols = 0
    found_types = 0
    found_long_names = 0
    for gid in gene_names:
        query_res: Dict[str, str] = gene_info_by_input.get(gid, dict())
        gene_symbol = query_res.get("symbol", "NA")
        if gene_symbol != "NA":
            found_symbols += 1
        gene_symbols.append(gene_symbol)
        gene_type = query_res.get("type_of_gene", "NA")
        if gene_type != "NA":
            found_types += 1
        gene_types.append(gene_type)
        gene_long_name = query_res.get("name", "NA")
        if gene_long_name != "NA":
            found_long_names += 1
        gene_long_names.append(gene_long_name)

    if found_symbols > 0.5 * len(gene_names):
        message("info", {"title": "Success", "body": f"Found {found_symbols} symbols"})
        print("Successfully mined for gene symbols")
    else:
        message(
            "warning",
            {
                "title": "Less than 50 Percent Gene Symbol Count",
                "body": f"Found {found_symbols} symbols",
            },
        )
        print("Low gene symbol count")
    if found_types > 0.5 * len(gene_names):
        message("info", {"title": "Success", "body": f"Found {found_types} types"})
        print("Successfully mined for gene types")
    else:
        message(
            "warning",
            {
                "title": "Less than 50 Percent Gene Type Count",
                "body": f"Found {found_types} types",
            },
        )
        print("Low gene type count")
    if found_long_names > 0.5 * len(gene_names):
        message(
            "info", {"title": "Success", "body": f"Found {found_long_names} long names"}
        )
        print("Successfully mined for gene long names")
    else:
        message(
            "warning",
            {
                "title": "Less than 50 Percent Gene Name Count",
                "body": f"Found {found_long_names} names",
            },
        )
        print("Low gene name count")

    message("info", {"title": f"Success", "body": ""})

    print("Generating h5ad files...")
    sample_h5ad_files: List[LatchFile] = []
    individual_h5ads = []
    for sample_name, quant_dir in zip(sample_names, quant_dirs):
        message(
            "info", {"title": f"Generating H5AD File For {sample_name}", "body": ""}
        )
        h5ad: anndata.AnnData = pyroe.load_fry(
            frydir=str(Path(quant_dir)), output_format="scRNA"
        )

        h5ad.obs_names = [x + f"_{sample_name}" for x in h5ad.obs_names]
        sample_annotations = [sample_name for x in h5ad.obs_names]
        h5ad.obs["sample"] = sample_annotations
        h5ad.var["mygene_symbol"] = gene_symbols
        h5ad.var["mygene_type"] = gene_types
        h5ad.var["mygene_name"] = gene_long_names

        h5ad.var["mygene_symbol"] = h5ad.var["mygene_symbol"].astype(str)
        h5ad.var["mygene_type"] = h5ad.var["mygene_type"].astype(str)
        h5ad.var["mygene_name"] = h5ad.var["mygene_name"].astype(str)

        cell_types = infer_cell_type(h5ad, mouse_tissue, human_tissue)
        if cell_types is not None:
            h5ad.obs["deepsort_celltype"] = [x[0] for x in cell_types]
            h5ad.obs["deepsort_subcelltype"] = [x[1] for x in cell_types]

        h5ad.write(f"{sample_name}_counts.h5ad")

        if remove_ambient_rna:
            run_decontX(Path(f"/root/{sample_name}_counts.h5ad"), sample_name)

        sample_h5ad_files.append(
            LatchFile(
                f"/root/{sample_name}_counts.h5ad",
                f"{output_name}{sample_name}/counts.h5ad",
            )
        )
        individual_h5ads.append(h5ad)

    combined_adata = anndata.concat(individual_h5ads)
    combined_adata.var["mygene_symbol"] = gene_symbols
    combined_adata.var["mygene_type"] = gene_types
    combined_adata.var["mygene_name"] = gene_long_names
    combined_adata.var["mygene_symbol"] = combined_adata.var["mygene_symbol"].astype(
        str
    )
    combined_adata.var["mygene_type"] = combined_adata.var["mygene_type"].astype(str)
    combined_adata.var["mygene_name"] = combined_adata.var["mygene_name"].astype(str)
    combined_adata.write(f"counts.h5ad")
    if remove_ambient_rna:
        run_decontX(Path("/root/counts.h5ad"), sample_name)

    message("info", {"title": f"Success", "body": ""})

    return (
        sample_h5ad_files,
        LatchFile(f"/root/counts.h5ad", f"{output_name}combined_counts.h5ad"),
    )


@large_task
def generate_report(
    map_dirs: List[LatchDir],
    permit_dirs: List[LatchDir],
    quant_dirs: List[LatchDir],
    sample_names: List[str],
    output_name: str,
) -> List[LatchFile]:
    reports: List[LatchFile] = []
    for sample_name, map_dir, permit_dir, quant_dir in zip(
        sample_names, map_dirs, permit_dirs, quant_dirs
    ):
        message(
            "info", {"title": f"Generating QC Report For {sample_name}", "body": ""}
        )
        print(f"Generating QC Report For {sample_name}...")
        report_cmd = [
            "Rscript",
            "qc.R",
            str(Path(map_dir)),
            str(Path(permit_dir)),
            str(Path(quant_dir)),
            quote(sample_name),
        ]

        report_process = subprocess.Popen(
            " ".join(report_cmd),
            shell=True,
            errors="replace",
            encoding="utf-8",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        while True:
            realtime_output = report_process.stdout.readline()

            if realtime_output == "" and report_process.poll() is not None:
                break
            if realtime_output:
                print("\t" + realtime_output.strip())

        retval = report_process.poll()
        if retval != 0:
            message(
                "error",
                {"title": "Alevin-Fry QC Report", "body": f"View logs to see error"},
            )
            raise ReportGenerationError(
                f"\tAlevin-Fry QC Report failed with exit code {retval}"
            )

        reports.append(
            LatchFile(
                f"/root/{sample_name}_alevinQC.html",
                f"{output_name}{sample_name}/alevinQC.html",
            )
        )

    message("info", {"title": "Success", "body": ""})
    print("QC Reports complete. Packaging Files...")
    return reports


@workflow
def scrnaseq(
    samples: List[Sample],
    read_length: Optional[int],
    run_name: str,
    latch_genome: LatchGenome,
    ta_ref_genome_fork: str = "database",
    celltype_fork: str = "default",
    output_location_fork: str = "default",
    technology: Technology = Technology.chromiumv3,
    custom_gtf: Optional[LatchFile] = None,
    custom_ref_genome: Optional[LatchFile] = None,
    splici_index: Optional[LatchDir] = None,
    custom_output_dir: Optional[LatchDir] = None,
    human_tissue: Optional[HumanTissue] = None,
    mouse_tissue: Optional[MouseTissue] = None,
    remove_ambient_rna: bool = False,
) -> Tuple[
    List[LatchDir],
    List[LatchDir],
    LatchFile,
    List[LatchFile],
    List[LatchFile],
]:
    """Performs alignment & quantification on Single Cell RNA-Sequencing reads.

    Single Cell RNA-Seq (Alignment & Quantification)
    ----

    This workflow will produce gene and transcript counts from single cell RNA-seq
    sample reads.

    This current iteration of the workflow has three steps:

    1. Generating an intron aware index [Splici Index](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2.supplementary-material)
    2. Selective alignment of reads using [Salmon](https://github.com/COMBINE-lab/salmon)
    3. Quantification of aligned reads using [Alevin-Fry](https://github.com/COMBINE-lab/alevin-fry)
    4. Quality Control using [alevinQC](https://github.com/csoneson/alevinQC)

    In each of these steps, the user is limited in the set of options available. In a few cases, such as the read length parameter
    for splici index creation, we estimate read length given the reads. In other cases, some functionality of the underlying tool
    is not exposed such as read geometry for the salmon alevin subcommand. If you need this functionality, please reach out to
    aidan@latch.bio so we can expose it for you.

    ### Output Format Specification

    U, S, and A refer to unspliced, spliced, ambiguous respectively.

    The output of this workflow is a directory containing the following files and subdirectories:
    * alevinQC.html: an alignment and quantification quality control report.
    * counts.h5ad: a gene by cell H5AD file containing summed spliced and ambiguous counts for each gene.
    * counts_velocity.h5ad: a gene by cell H5AD file containing the "spliced" layer, which contains the S+A counts, and the "unspliced" layer, which contains the U counts.
    * counts_USA.h5ad: a gene by cell H5AD file containing a layer for each count type.
    * raw_counts/: a directory containing the raw output from alevin-fry quantification.
    * intermediate_mapping/: a directory containing the intermediate mapping RAD files from salmon alevin. Refer to [Splici Index](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2.supplementary-material) for further information.
    * quant_preprocess/: a directory containing the collated RAD files along with the binary cellular barcode whitelist files.
    * splici_index/: a directory containing the splici index files. This will only be generated if a splici index is not provided.

    ### Authors and Maintainers of Underlying Tools

    Credit is due for the splici index, selctive alignment (Salmon), and quantification (Alevin-Fry) to the [COMBINELab](https://combine-lab.github.io)
    who have created a toolset for the analysis of single cell and bulk RNA-seq data. These are the tools used in this workflow and all intellectual
    credit is theirs.

    Credit is also due to the author of [alevinQC](https://github.com/csoneson/alevinQC) for collecting the data from alignment and quantification
    into a single HTML report.

    ### More about Selective Alignment

    This workflow uses Salmon's selective alignment described in this
    [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8),
    which achieves greater accuracy than traditional alignment methods while
    using fewer computational resources.

    ### Justification for the Splici Index

    "The term splici is shorthand for spliced + intronic reference sequence. This reference sequence is prepared
    by extracting the spliced transcripts from the reference genome according to the desired annotation (e.g. unfiltered, or filtered for certain
    classes of transcripts), as well as the collapsed intervals corresponding to the introns of genes. The intronic sequences of the splici reference
    play important roles in the various kinds of experiments discussed in this paper. For scRNA-seq data, although one typically focuses on the fragments
    and UMIs arising from the (spliced) transcriptome,and only considers the spliced and ambiguous counts when performing downstream analyses, the
    intronic sequences act similarly to decoy sequences proposed by Srivastava et al. (4). They account for fragments that might otherwise selectively-align
    or pseudoalign to the transcriptome with lower quality, but that, in fact, derive from some unspliced RNA transcript molecule."
    -- [From the supplementary material section of the alevin fry paper.](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v2.supplementary-material

    __metadata__:
        display_name: Single Cell RNA-seq
        author:
            name: LatchBio
            email: aidan@latch.bio
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
                  using less computational resources.

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

        - section: Automatic Celltype Annotation
          flow:
          - text: >-
                  This workflow uses [scDeepSort](https://www.biorxiv.org/content/10.1101/2020.05.13.094953v1) to automatically add celltype annotations
                  to cells from select human and mouse tissues. In the future, we plan on supporting additional organisms and tissues. The scDeepSort algorithm
                  is trained on Tabula Muris Senis (TMS) and Tabula Sapiens (TS) data. The list of supported tissues is available [here](https://github.com/ZJUFanLab/scDeepSort/wiki/Mouse-tissues-and-cell-types).
                  Enabling this option will add two pieces of metadata to each cell: `deepsort_celltype` and `deepsort_subcelltype`.
          - fork: celltype_fork
            flows:
                default:
                    display_name: None
                    flow:
                    - text:
                        No celltype annotation will be performed.
                human:
                    display_name: Human
                    _tmp_unwrap_optionals:
                        - human_tissue
                    flow:
                    - params:
                        - human_tissue
                mouse:
                    display_name: Mouse
                    _tmp_unwrap_optionals:
                        - mouse_tissue
                    flow:
                    - params:
                        - mouse_tissue

        - section: Ambient RNA Removal
          flow:
            - text: >-
                  This workflow uses [decontX](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6)
                  to remove the influence of contaminating ambient RNA from the
                  counts matrix. In short, the algorithm explicitly defines the contamination distribution to be a weighted
                  combination of all other cell population distributions, allowing it to
                  tease out what ambient RNAs are floating around in the sample and correct all cells
                  counts for them. Enabling this parameter will modify the counts and add two
                  pieces of metadata to each cell: `decontX_contamination` and `decontX_clusters`.
            - params:
                - remove_ambient_rna

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
                        viewer - SCRNA-Seq Outputs/"Run Name"
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

        alignment_quantification_tools:

          __metadata__:
            display_name: Alignment & Quantification Method

        ta_ref_genome_fork:

            __metadata__:
                display_name: Reference Genome Source

        celltype_fork:

        human_tissue:
            __metadata__:
                display_name: Human Tissue

        mouse_tissue:
            __metadata__:
                display_name: Mouse Tissue

        mouse_tissue:

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
                detail: (gtf / gff, text or gzipped))

        splici_index:
            You are able to provide a prebuilt splici index directory for your genome.
            This will speed up run time as the index is generated if none is provided.
            This should be the output from a previous scrnaseq run or generated using
            this [tutorial](https://combine-lab.github.io/alevin-fry-tutorials/2021/improving-txome-specificity/).

          __metadata__:
            display_name: Splici index

        read_length:
            The length of the reads in the sample. This is used to generate the splici index.
            If you are providing an index, you can ignore this parameter. If you do not provide
            an index or fill out this parameter, the read length will be estimated from the
            first 100 reads in the first sample.
            __metadata__:
                display_name: Reads Length
                batch_table_column: true

        save_indices:
            If you provided a custom genome you can output the alignment
            indexes generated from this run for use in future runs. This will
            speed up runtime since the workflow doesn't have to then regenerate
            the indexes.

          __metadata__:
            display_name: Save Generated Reference Indexes

        remove_ambient_rna:
            Enabling this parameter will modify the counts and add two
            pieces of metadata to each cell: `decontX_contamination` and `decontX_clusters`.

          __metadata__:
            display_name: Run DecontX

        run_name:
          A name for this analysis run, this will be used to name outputs from
          this run.

          __metadata__:
            batch_table_column: true
            display_name: Run Name

        output_location_fork:
            __metadata__:
                display_name: Output Location

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

    (mapping_dirs, sample_names) = map_reads(
        splici_index=splici_index,
        samples=samples,
        technology=technology,
        output_name=output_name,
    )

    (preprocessed_quant_dirs, quantified_read_dirs) = quantify_reads(
        mapping_dirs=mapping_dirs,
        sample_names=sample_names,
        tg_map=t2g,
        output_name=output_name,
    )

    (individual_counts, combined_counts) = h5ad(
        quant_dirs=quantified_read_dirs,
        sample_names=sample_names,
        output_name=output_name,
        human_tissue=human_tissue,
        mouse_tissue=mouse_tissue,
        remove_ambient_rna=remove_ambient_rna,
    )

    reports = generate_report(
        map_dirs=mapping_dirs,
        permit_dirs=preprocessed_quant_dirs,
        quant_dirs=quantified_read_dirs,
        sample_names=sample_names,
        output_name=output_name,
    )

    return (
        preprocessed_quant_dirs,
        quantified_read_dirs,
        combined_counts,
        individual_counts,
        reports,
    )


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
