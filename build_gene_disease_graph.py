import pandas as pd
import networkx as nx
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
REL_PATH = BASE_DIR / "relationships.tsv"
HGNC_PATH = BASE_DIR / "hgnc_complete_set.tsv"

def main():
    df = pd.read_csv(REL_PATH, sep="\t")

    gd1 = df[
        (df["Entity1_type"] == "Gene") &
        (df["Entity2_type"] == "Disease")
    ][[
        "Entity1_id", "Entity1_name",
        "Entity2_id", "Entity2_name",
        "Evidence", "Association", "PK", "PD", "PMIDs"
    ]].copy()

    gd1.columns = [
        "gene_id", "gene_symbol",
        "disease_id", "disease_name",
        "evidence_type", "association", "pk_related", "pd_related", "pmids"
    ]

    gd2 = df[
        (df["Entity1_type"] == "Disease") &
        (df["Entity2_type"] == "Gene")
    ][[
        "Entity2_id", "Entity2_name",
        "Entity1_id", "Entity1_name",
        "Evidence", "Association", "PK", "PD", "PMIDs"
    ]].copy()

    gd2.columns = [
        "gene_id", "gene_symbol",
        "disease_id", "disease_name",
        "evidence_type", "association", "pk_related", "pd_related", "pmids"
    ]

    edges_df = pd.concat([gd1, gd2], ignore_index=True)

    # Remove Duplicate
    edges_df["pair_key"] = edges_df.apply(
        lambda r: tuple(sorted([str(r["gene_id"]), str(r["disease_id"])])),
        axis=1
    )
    edges_df = edges_df.drop_duplicates(subset=["pair_key", "association", "evidence_type"])
    edges_df = edges_df.drop(columns=["pair_key"])

    # =========================
    # Add gene family from HGNC
    # =========================
    hgnc_df = pd.read_csv(HGNC_PATH, sep="\t", dtype=str)

    gene_family_df = hgnc_df[["symbol", "gene_group"]].copy()
    gene_family_df.columns = ["gene_symbol", "gene_family"]

    # normalize key for safer matching
    gene_family_df["merge_symbol"] = gene_family_df["gene_symbol"].astype(str).str.strip().str.upper()
    edges_df["merge_symbol"] = edges_df["gene_symbol"].astype(str).str.strip().str.upper()

    # remove duplicated HGNC symbols if any
    gene_family_df = gene_family_df.drop_duplicates(subset=["merge_symbol"])

    edges_df = edges_df.merge(
        gene_family_df[["merge_symbol", "gene_family"]],
        on="merge_symbol",
        how="left"
    )

    edges_df = edges_df.drop(columns=["merge_symbol"])

    # Fill missing family
    edges_df["gene_family"] = edges_df["gene_family"].fillna("Unknown")

    print("Number of gene-disease edges:", len(edges_df))
    print(edges_df.head())

    edges_df.to_csv("gene_disease_edges.csv", index=False)
    print("Saved to gene_disease_edges.csv")

    # Create Graph
    G = nx.Graph()

    for _, row in edges_df.iterrows():
        G.add_node(
            row["gene_symbol"],
            node_type="gene",
            pharmgkb_id=row["gene_id"],
            gene_family=row["gene_family"]
        )
        G.add_node(
            row["disease_name"],
            node_type="disease",
            pharmgkb_id=row["disease_id"]
        )

        G.add_edge(
            row["gene_symbol"],
            row["disease_name"],
            evidence_type=row["evidence_type"],
            association=row["association"],
            pk_related=row["pk_related"],
            pd_related=row["pd_related"],
            pmids=row["pmids"],
            gene_family=row["gene_family"]
        )

    print("Number of nodes:", G.number_of_nodes())
    print("Number of edges in graph:", G.number_of_edges())

    gene_degrees = []
    for node, degree in G.degree():
        if G.nodes[node].get("node_type") == "gene":
            gene_degrees.append((node, degree))

    gene_degrees = sorted(gene_degrees, key=lambda x: x[1], reverse=True)
    print("\nTop 10 genes with most disease links:")
    for gene, deg in gene_degrees[:10]:
        print(f"{gene}: {deg}")

if __name__ == "__main__":
    main()