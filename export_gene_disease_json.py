import json
from pathlib import Path
import pandas as pd

# ------------------------------------------------------------
# Export gene-disease CSV to knowledge graph JSON
# Output format:
# {
#   "metadata": {...},
#   "nodes": [...],
#   "links": [...]
# }
# ------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
INPUT_CSV = BASE_DIR / "gene_disease_edges.csv"
OUTPUT_JSON = BASE_DIR / "gene_disease.json"


def clean_value(value):
    """Convert pandas/CSV missing values to None and strip strings."""
    if pd.isna(value):
        return None
    if isinstance(value, str):
        value = value.strip()
        return value if value else None
    return value


def parse_pmids(value):
    """Turn PMID text into a list of strings."""
    value = clean_value(value)
    if value is None:
        return []

    text = str(value)
    for sep in [";", "|"]:
        text = text.replace(sep, ",")

    pmids = [item.strip() for item in text.split(",") if item.strip()]
    return pmids


def build_node_id(node_type, pharmgkb_id, label):
    """Create a stable node id for graph links."""
    if pharmgkb_id:
        return f"{node_type}:{pharmgkb_id}"
    safe_label = str(label).strip().replace(" ", "_")
    return f"{node_type}:{safe_label}"


def get_family_value(row):
    """Read family from either 'family' or 'gene_family' column."""
    if "family" in row.index:
        value = clean_value(row["family"])
        if value is not None:
            return value

    if "gene_family" in row.index:
        value = clean_value(row["gene_family"])
        if value is not None:
            return value

    return None


def main():
    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)

    required_columns = {
        "gene_id", "gene_symbol", "disease_id", "disease_name",
        "evidence_type", "association", "pk_related", "pd_related", "pmids"
    }
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    nodes = {}
    links = []

    for _, row in df.iterrows():
        gene_id = clean_value(row["gene_id"])
        gene_symbol = clean_value(row["gene_symbol"])
        disease_id = clean_value(row["disease_id"])
        disease_name = clean_value(row["disease_name"])
        family = get_family_value(row)

        if gene_symbol is None or disease_name is None:
            continue

        gene_node_id = build_node_id("gene", gene_id, gene_symbol)
        disease_node_id = build_node_id("disease", disease_id, disease_name)

        # Add gene node
        if gene_node_id not in nodes:
            nodes[gene_node_id] = {
                "id": gene_node_id,
                "label": gene_symbol,
                "type": "gene",
                "group": "gene",
                "pharmgkb_id": gene_id,
                "family": family,
            }

        # Add disease node
        if disease_node_id not in nodes:
            nodes[disease_node_id] = {
                "id": disease_node_id,
                "label": disease_name,
                "type": "disease",
                "group": "disease",
                "pharmgkb_id": disease_id,
            }

        # Add link
        links.append({
            "source": gene_node_id,
            "target": disease_node_id,
            "relation": "gene-disease",
            "evidence_type": clean_value(row["evidence_type"]),
            "association": clean_value(row["association"]),
            "pk_related": clean_value(row["pk_related"]),
            "pd_related": clean_value(row["pd_related"]),
            "pmids": parse_pmids(row["pmids"]),
            "family": family,
        })

    graph_json = {
        "metadata": {
            "graph_type": "gene-disease knowledge graph",
            "input_file": INPUT_CSV.name,
            "node_count": len(nodes),
            "link_count": len(links),
        },
        "nodes": list(nodes.values()),
        "links": links,
    }

    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(graph_json, f, indent=2, ensure_ascii=False)

    print(f"Saved JSON to: {OUTPUT_JSON}")
    print(f"Nodes: {len(nodes)}")
    print(f"Links: {len(links)}")


if __name__ == "__main__":
    main()