import csv
import os
import re
import shutil
import tempfile
from collections import defaultdict
from Bio import SeqIO


def select_best_isoforms(pp_file, score_cols=None):
    """
    Select the best protein variants from a PhyloProfile (pp) file.

    Selection rules
    ---------------
    - Variants labelled exactly 's' are always retained.
    - Self isoforms ('si', 'sio', ...) compete with each other.
        * If the best self isoform has a higher score than 's',
          both 's' and the best self isoform(s) are retained.
        * Otherwise only 's' is retained.
    - Numbered orthologs ('1', '1i', '1s', '2', '2so', ...)
      compete within their ortholog number.
      Only the highest scoring variant(s) are retained.
    - Ties are preserved.

    Parameters
    ----------
    pp_file : str
        Input PhyloProfile file.

    score_cols : list[str], optional
        Score columns used for ranking.
        If None, all columns except geneID, ncbiID and orthoID
        are averaged.

    Returns
    -------
    set
        Set of orthoID entries to retain.
    """

    groups = defaultdict(lambda: {
        "self": None,
        "self_iso": None,
        "best": {}
    })

    with open(pp_file) as f:
        reader = csv.DictReader(
            f,
            delimiter="\t"
        )
        if score_cols is None:
            excluded = {
                "geneID",
                "ncbiID",
                "orthoID"
            }
            score_cols = [
                c for c in reader.fieldnames
                if c not in excluded
            ]

        for row in reader:
            fields = row["orthoID"].split("|")
            if len(fields) != 4:
                continue
            gene_id = row["geneID"]
            group = fields[1]
            # combine geneID and ortholog group
            group_key = (gene_id, group)

            protein = fields[2]
            orthotype = fields[3]
            score = sum(
                float(row[c])
                for c in score_cols
            ) / len(score_cols)

            # reference sequence
            if orthotype == "s":
                groups[group_key]["self"] = (
                    row["orthoID"],
                    protein,
                    score
                )
                continue

            # self isoforms
            if orthotype.startswith("s"):
                current = groups[group_key]["self_iso"]
                if current is None:
                    groups[group_key]["self_iso"] = (
                        score,
                        [(row["orthoID"], protein)]
                    )
                else:
                    best_score, proteins = current
                    if score > best_score:
                        groups[group_key]["self_iso"] = (
                            score,
                            [(row["orthoID"], protein)]
                        )
                    elif score == best_score:
                        proteins.append(
                            (row["orthoID"], protein)
                        )
                continue

            # numbered orthologs
            match = re.match(
                r"\d+",
                orthotype
            )

            if match is None:
                continue
            num = match.group()
            if num not in groups[group_key]["best"]:
                groups[group_key]["best"][num] = (
                    score,
                    [(row["orthoID"], protein)]
                )
            else:
                best_score, proteins = groups[group_key]["best"][num]
                if score > best_score:
                    groups[group_key]["best"][num] = (
                        score,
                        [(row["orthoID"], protein)]
                    )
                elif score == best_score:
                    proteins.append(
                        (row["orthoID"], protein)
                    )

    keep_ids = set()
    for (gene_id, group), data in groups.items():
        # always keep s
        if data["self"] is not None:
            self_id, _, self_score = data["self"]
            keep_ids.add(self_id)
            # keep best si/sio only if better than s
            if data["self_iso"] is not None:
                iso_score, isoforms = data["self_iso"]
                if iso_score > self_score:
                    for ortho_id, _ in isoforms:
                        keep_ids.add(ortho_id)

        # numbered orthologs
        for _, (_, proteins) in data["best"].items():
            for ortho_id, _ in proteins:
                keep_ids.add(ortho_id)
    return keep_ids


def temporary_output(input_file, output_file):
    """
    Create a temporary file if overwriting is required.
    """
    if os.path.abspath(input_file) == os.path.abspath(output_file):
        tmp = tempfile.NamedTemporaryFile(
            delete=False,
            dir=os.path.dirname(input_file)
        )
        tmp.close()
        return tmp.name, True
    return output_file, False


def finalize_output(tmp_file, output_file, replace):
    """
    Replace original file after successful writing.
    """
    if replace:
        shutil.move(
            tmp_file,
            output_file
        )


def filter_fasta(input_file, output_file, keep_ids):
    """
    Keep only selected FASTA sequences.
    Output sequences are written in one line,
    preserving the input FASTA style.
    """
    with open(output_file, "w") as out:
        for record in SeqIO.parse(
            input_file,
            "fasta"
        ):
            if record.id in keep_ids:
                out.write(
                    f">{record.id}\n"
                )
                out.write(
                    str(record.seq) + "\n"
                )


def filter_domain_file(input_file, output_file, keep_ids):
    """
    Filter domain annotation files by orthoID.
    """
    with open(input_file) as inp, open(output_file, "w") as out:
        for line in inp:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            if fields[1] in keep_ids:
                out.write(line)


def filter_pp_file(input_file, output_file, keep_ids):
    """
    Filter PhyloProfile file by orthoID.
    """
    with open(input_file) as inp, open(output_file, "w") as out:
        header = inp.readline()
        out.write(header)
        for line in inp:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            if fields[2] in keep_ids:
                out.write(line)


def get_output_path(input_file, output_dir):
    """
    Return output path.
    """
    if output_dir is None:
        return input_file
    return os.path.join(
        output_dir,
        os.path.basename(input_file)
    )


def filter_isoforms(
        pp_file,
        fasta_file = None,
        forward_file = None,
        reverse_file = None,
        output_dir = None,
        score_cols = None):
    """
    Select best isoforms and filter all related files.
    """

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)


    keep_ids = select_best_isoforms(
        pp_file,
        score_cols=score_cols
    )


    files = [
        (pp_file, filter_pp_file),
        (fasta_file, filter_fasta),
        (forward_file, filter_domain_file),
        (reverse_file, filter_domain_file)
    ]


    for input_file, func in files:
        # Skip optional files that are not provided or missing
        if input_file is None:
            continue

        if not os.path.isfile(input_file):
            print(
                f"Warning: file not found, skipping: {input_file}"
            )
            continue

        output_file = get_output_path(
            input_file,
            output_dir
        )

        tmp_file, replace = temporary_output(
            input_file,
            output_file
        )

        func(
            input_file,
            tmp_file,
            keep_ids
        )

        finalize_output(
            tmp_file,
            output_file,
            replace
        )

    return keep_ids