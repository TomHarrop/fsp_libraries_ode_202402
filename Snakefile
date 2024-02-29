#!/usr/bin/env python3

from functools import cache
from pathlib import Path
from snakemake.logging import logger
import pandas as pd
import pickle
import re
import tempfile

#############
# FUNCTIONS #
#############


@cache
def check_pool_file(pool_file):
    return Path(pool_file).is_file()


@cache
def find_sample_input(wildcards):
    my_plate = sample_to_plate[wildcards.sample]
    return Path(
        outdir,
        "010_tcdemux_unpooled",
        my_plate,
        f"{wildcards.sample}.{{read}}.fastq.gz",
    )


@cache
def get_all_pools(sample_data_file):
    return sorted(set(read_sample_csv(sample_data_file)["pool_name"]))


@cache
def get_all_samples(sample_data_file):
    return sorted(set(read_sample_csv(sample_data_file).index))


@cache
def read_sample_csv(sample_data_file):
    my_csv = pd.read_csv(sample_data_file, index_col="name")
    for i, row in my_csv.iterrows():
        sanitised_name = sanitise_sample_name(i)
        if sanitised_name != i:
            logger.warning(f"Replacing sample name {i} with {sanitised_name}")
            my_csv.rename(index={i: sanitised_name}, inplace=True)
    return my_csv


@cache
def sanitise_sample_name(sample_name):
    return re.sub("[^a-zA-Z0-9_\-]", "", sample_name)


###########
# GLOBALS #
###########

raw_read_directory = Path("data", "raw_reads")
outdir = Path("output")
logdir = Path(outdir, "logs")

# set up a temporary directory for this run
try:
    run_tmpdir = config["run_tmpdir"]
    print(f"Caught run_tmpdir {run_tmpdir}")
except KeyError as e:
    print(f"{e} not set in config")
    run_tmpdir = tempfile.mkdtemp()
    print(f"Setting run_tmpdir to {run_tmpdir}")
    print("This probably won't work on a cluster!")


adaptor_files = [Path("data", "adaptors", "bbmap_39.01_adaptors.fa")]

# containers
tcdemux = "docker://quay.io/biocontainers/tcdemux:0.0.24--pyhdfd78af_0"

# modules
module_tag = "0.0.42"
hybpiper_snakefile = github(
    "tomharrop/smk-modules",
    path="modules/hybpiper/Snakefile",
    tag=module_tag,
)

########
# MAIN #
########

all_plates = [
    Path(x) for x in raw_read_directory.glob("*") if Path(x).is_dir()
]

sample_dict = {}
pool_dict = {}
sample_to_plate = {}

for plate_path in all_plates:
    pool_dict[plate_path.name] = get_all_pools(
        Path(plate_path, "sample_data.csv")
    )
    my_samples = get_all_samples(Path(plate_path, "sample_data.csv"))
    sample_dict[plate_path.name] = my_samples
    for sample in my_samples:
        sample_to_plate[sample] = plate_path.name

all_samples = sorted(set(sample_to_plate.keys()))


# This is necessary because some samples have zero reads. The pipeline has to
# be run in two steps for now.
@cache
def sample_status_file_exists(status_file):
    return Path(status_file).is_file()


def check_sample_status(status_file):
    status_path = Path(status_file).resolve()
    with open(status_path, "rb") as f:
        return pickle.load(f)


def get_hybpiper_input(all_samples):
    hybpiper_samples = []
    for sample in all_samples:
        status_file = Path(
            outdir,
            "010_tcdemux_unpooled",
            "all_samples",
            f"{sample}.check",
        ).resolve()
        if sample_status_file_exists(status_file):
            if check_sample_status(sample):
                hybpiper_samples.append(sample)
            else:
                logger.warning(f"Dropping sample {sample}")
    return sorted(set(hybpiper_samples))


hybpiper_samples = get_hybpiper_input(all_samples)

# see tomharrop/gap_library_tests/data/target_files/README.md for how this was
# generated
target_file = Path("data", "target_files", "gap_targets.fa.gz")


#########
# RULES #
#########


# rule target:
#     input:
#         expand(
#             Path(
#                 outdir,
#                 "010_tcdemux_unpooled",
#                 "all_samples",
#                 "{sample}.{read}.fastq.gz",
#             ),
#             sample=all_samples,
#             read=["r1", "r2"],
#         ),


if len(hybpiper_samples) > 0:

    module hybpiper:
        snakefile:
            hybpiper_snakefile
        config:
            {
                "outdir": Path(outdir, "020_hybpiper"),
                "read_directory": Path(
                outdir, "010_tcdemux_unpooled", "all_samples"
                ),
                "run_tmpdir": run_tmpdir,
                "sample_list": get_hybpiper_input(all_samples),
                "target_file": target_file,
            }

    use rule * from hybpiper as hybpiper_*

    use rule hybpiper_assemble from hybpiper as hybpiper_hybpiper_assemble with:
        resources:
            time=lambda wildcards, attempt: 240 * attempt,
            mem_mb=lambda wildcards, attempt: 16e3 * attempt,

else:
    logger.warning(
        """No demuxed samples detected.
Hybpiper target will be skipped.
This warning should go away after you run demux_target."""
    )


rule demux_target:
    input:
        expand(
            Path(
                outdir,
                "010_tcdemux_unpooled",
                "all_samples",
                "{sample}.check",
            ),
            sample=all_samples,
        ),


rule check_demuxed_samples:
    input:
        r1=Path(
            outdir,
            "010_tcdemux_unpooled",
            "all_samples",
            "{sample}.r1.fastq.gz",
        ),
        r2=Path(
            outdir,
            "010_tcdemux_unpooled",
            "all_samples",
            "{sample}.r2.fastq.gz",
        ),
    output:
        Path(
            outdir,
            "010_tcdemux_unpooled",
            "all_samples",
            "{sample}.check",
        ),
    threads: 1
    resources:
        mem_mb=int(1e3),
        time=1,
    run:
        sample_status = (Path(input.r1).stat().st_size > 20) and (
            Path(input.r2).stat().st_size > 20
        )
        with open(output[0], "wb") as f:
            pickle.dump(sample_status, f)



rule collect_demuxed_files:
    input:
        find_sample_input,
    output:
        Path(
            outdir,
            "010_tcdemux_unpooled",
            "all_samples",
            "{sample}.{read}.fastq.gz",
        ),
    threads: 1
    resources:
        mem_mb=int(1e3),
        time=1,
    shell:
        "ln -s "
        "$( readlink -f {input} ) "
        "$( readlink -f {output} )"


for plate_path in all_plates:
    my_plate = plate_path.name
    my_adaptors = adaptor_files + [Path(plate_path, "adapters.fa")]

    rule:
        input:
            read_files=[x for x in Path(plate_path).glob("*.fastq*")],
            sample_data=Path(run_tmpdir, my_plate, "sample_data.csv"),
            adaptor_files=my_adaptors,
        output:
            expand(
                Path(
                    outdir,
                    "010_tcdemux_unpooled",
                    my_plate,
                    "{sample}.{read}.fastq.gz",
                ),
                sample=sample_dict[my_plate],
                read=["r1", "r2", "unpaired"],
            ),
        params:
            outdir=lambda wildcards, output: Path(output[0]).resolve().parent,
            read_dir=lambda wildcards, input: Path(input.read_files[0])
            .resolve()
            .parent,
            restart_times=5,
            mem_gb=lambda wildcards, resources: int(
                resources.mem_mb / 1e3 * 0.9
            ),
        log:
            Path(logdir, f"tcdemux_unpooled.{my_plate}.log"),
        resources:
            mem_mb=lambda wildcards, threads, attempt: (threads // 5)
            * 8e3
            * attempt,
            time=lambda wildcards, attempt: 2880 * attempt,
        threads: 60
        container:
            tcdemux
        shell:
            "tcdemux "
            "--sample_data {input.sample_data} "
            "--adaptors {input.adaptor_files} "
            "--outdir {params.outdir} "
            "--read_directory {params.read_dir} "
            "--threads {threads} "
            "--mem_gb {params.mem_gb} "
            "--restart_times {params.restart_times} "
            "&> {log}"

    rule:
        input:
            Path(plate_path, "sample_data.csv"),
        output:
            Path(run_tmpdir, my_plate, "sample_data.csv"),
        params:
            plate_path,
        run:
            sample_csv = read_sample_csv(input[0])
            for i, row in sample_csv.iterrows():
                if not (
                    check_pool_file(Path(params[0], row.r1_file))
                    and check_pool_file(Path(params[0], row.r2_file))
                ):
                    raise ValueError(f"Missing pool file for {i}")
            sample_csv.to_csv(output[0])
