from collections import defaultdict
import re


def parse_gff_protein_to_gene(gff_file):
    """
    Returns:
        protein_to_gene: dict
            protein_id -> GeneID
    """
    protein_to_gene = {}

    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue

            feature_type = fields[2]
            if feature_type != "CDS":
                continue

            attrs = fields[8]

            protein_match = re.search(r"protein_id=([^;]+)", attrs)
            gene_match = re.search(r"GeneID:(\d+)", attrs)

            if protein_match and gene_match:
                protein_id = protein_match.group(1)
                gene_id = gene_match.group(1)

                protein_to_gene[protein_id] = gene_id

    return protein_to_gene


def find_isoforms(protein_to_gene, protein_id):
    """
    Return all isoforms (protein IDs) belonging to the same gene as protein_id.

    Parameters
    ----------
    protein_to_gene : dict
        Dictionary mapping protein IDs to gene IDs.
        Example: {'XP_040513265.1': '121108972', ...}

    protein_id : str
        Protein ID for which to find isoforms.

    Returns
    -------
    list
        List of protein IDs sharing the same gene ID.
    """

    # Check if protein exists
    if protein_id not in protein_to_gene:
        raise ValueError(f"{protein_id} not found in protein_to_gene dictionary")

    # Get the gene ID of the query protein
    gene_id = protein_to_gene[protein_id]

    # Find all proteins with the same gene ID
    isoforms = [
        prot for prot, gene in protein_to_gene.items()
        if gene == gene_id
    ]

    return isoforms


def build_protein_group_map(proteins, protein_to_gene):
    """
    Returns:
        protein_to_group : dict
            {protein_id: gene_id}
    """
    protein_to_group = {}

    for prot in proteins:
        gene_id = protein_to_gene.get(prot)

        if gene_id is None:
            gene_id = f"UNKNOWN:{prot}"

        protein_to_group[prot] = gene_id

    return protein_to_group


def assign_orthotypes(protein_ids, protein_to_gene, seed_id, filterIsoforms):
    """
    Parameters
    ----------
    protein_ids : list[str]
        Ordered list of proteins

    protein_to_gene : dict
        protein_id -> gene_id

    seed_id: str
        Seed protein ID: to specify orthologs that are indeed isoforms with seed 

    filterIsoforms: str
        all_isoforms: get all isoforms for each protein in protein_ids
        fas_isoform: get all isoforms for each protein in protein_ids (and filter later by best FAS score)

    Returns
    -------
    dict
        protein_id -> orthotype (e.g. "1", "1i", "2", "2i", ...)
    """

    gene_to_index = {}   # gene_id -> numeric index
    seen_gene = set()
    result = {}
    other_isoform = {}
    gene_counter = 0

    # get gene ID of seed
    seed_gene = ""
    if protein_to_gene:
        seed_gene = protein_to_gene.get(seed_id)

    for prot in protein_ids:
        gene = protein_to_gene.get(prot)
        # handle unknowns explicitly if needed
        if gene is None:
            gene = f"UNKNOWN:{prot}"
        # assign gene index if first time seen
        if gene not in gene_to_index:
            gene_counter += 1
            gene_to_index[gene] = gene_counter
            is_first_in_gene = True
            other_isoform = {}
        else:
            is_first_in_gene = False

        idx = gene_to_index[gene]
        if gene == seed_gene:
            idx = "s"
        # assign orthotype
        if is_first_in_gene:
            result[prot] = str(idx)
            if idx == "s":
                result[prot] = "si"
        else:
            result[prot] = f"{idx}i"

        if filterIsoforms == "all_isoforms" or filterIsoforms == "fas_isoform":
            if protein_to_gene:
                # get other isoforms for this protein
                isoforms = find_isoforms(protein_to_gene, prot)
                for isoform in isoforms:
                    if isoform not in protein_ids:
                        result[isoform] = f"{idx}io"
                        other_isoform[isoform] = gene
    return result
