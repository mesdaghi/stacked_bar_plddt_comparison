import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_and_prepare(old_file, new_file=None):
    # Case 1: single file with source column (_old / _new)
    if new_file is None:
        df = pd.read_csv(old_file, sep="\t")
        before = df[df["source"].str.endswith("_old")].copy()
        after = df[df["source"].str.endswith("_new")].copy()
    else:
        # Case 2: two separate files
        before = pd.read_csv(old_file, sep="\t")
        after = pd.read_csv(new_file, sep="\t")

    # Gene ID extraction
    def extract_gene_id(df, esm=False):
        if "gene_name" in df.columns:  # guard in case column name is different
            if esm:
                df = df[df["gene_name"] != "gene_name"]
                gene_id = df["gene_name"].str.extract(r'name=([^ ]+)')[0].fillna(
                    df["gene_name"].str.split().str[0]
                )
                gene_id = gene_id.str.rsplit("-", n=1).str[0]
            else:
                gene_id = df["gene_name"].str.split("-").str[0]
        elif "gene" in df.columns:
            gene_id = df["gene"].str.split("-").str[0]
        else:
            raise ValueError("Expected 'gene_name' or 'gene' column not found.")

        dup_ids = gene_id[gene_id.duplicated(keep=False)]
        gene_id[gene_id.isin(dup_ids)] = None
        df = df.assign(gene_id=gene_id).dropna(subset=["gene_id"])
        return df

    esm_flag = (
        (old_file and "esm" in old_file.lower())
        or (new_file and "esm" in new_file.lower())
    )
    before = extract_gene_id(before, esm=esm_flag)
    after = extract_gene_id(after, esm=esm_flag)

    print(f"{old_file} first gene IDs:", before["gene_id"].head(10).tolist())
    if new_file:
        print(f"{new_file} first gene IDs:", after["gene_id"].head(10).tolist())

    merged = pd.merge(before, after, on="gene_id", suffixes=("_old", "_new"))
    print(f"Loaded {old_file} → {merged.shape[0]} merged rows")

    # Convert cumulative PLDDT columns to numeric
    for prefix in ["old", "new"]:
        for t in [50, 60, 70, 80, 90]:
            merged[f"residue_plddt_{t}_{prefix}"] = pd.to_numeric(
                merged[f"residue_plddt_{t}_{prefix}"], errors="coerce"
            )

        # Compute ≥ thresholds using absolute counts (same as Protenix)
        merged[f"≥50_{prefix}"] = (
            merged[f"residue_plddt_50_{prefix}"]
            + merged[f"residue_plddt_60_{prefix}"]
            + merged[f"residue_plddt_70_{prefix}"]
            + merged[f"residue_plddt_80_{prefix}"]
            + merged[f"residue_plddt_90_{prefix}"]
        )
        merged[f"≥60_{prefix}"] = (
            merged[f"residue_plddt_60_{prefix}"]
            + merged[f"residue_plddt_70_{prefix}"]
            + merged[f"residue_plddt_80_{prefix}"]
            + merged[f"residue_plddt_90_{prefix}"]
        )
        merged[f"≥70_{prefix}"] = (
            merged[f"residue_plddt_70_{prefix}"]
            + merged[f"residue_plddt_80_{prefix}"]
            + merged[f"residue_plddt_90_{prefix}"]
        )
        merged[f"≥80_{prefix}"] = (
            merged[f"residue_plddt_80_{prefix}"]
            + merged[f"residue_plddt_90_{prefix}"]
        )
        merged[f"≥90_{prefix}"] = merged[f"residue_plddt_90_{prefix}"]

    return merged


# --- Load datasets ---

# Protenix
protenix_df = load_and_prepare(
    "../protenix/before2.tsv",
    "../protenix/after2.tsv"
)

# Protenix_1
protenix1_df = load_and_prepare(
    "../protenix/protenix_1/protenix_1_before.tsv",
    "../protenix/protenix_1/protenix_1_after.tsv"
)

# ESM
esm_df = load_and_prepare(
    "../fumigatus/ox_461_esm_before.tsv",  # before
    "../fumigatus/ox_461_esm_after.tsv"  # after
)

# AF3 (single file with `source` column)
af3_df = load_and_prepare("AfumAf293_AF3_long.tsv")


# Load 461 proteins
proteins_461 = pd.read_csv("../fumigatus/461.csv")['old_id'].astype(str).str.upper().tolist()
proteins_461_set = set(proteins_461)

# List of DataFrames and method names
dfs = [
#    (protenix_df, "Protenix"),
    (protenix1_df, "Protenix_1"),
    (esm_df, "ESM"),
    (af3_df, "AF3"),
]

summary = []

for df, method in dfs:
    # Determine protein ID column
    if 'gene_old' in df.columns:
        id_col = 'gene_old'
    elif 'gene_name_old' in df.columns:
        id_col = 'gene_name_old'
    else:
        id_col = df.columns[0]  # fallback

    # Normalize IDs: remove suffix after '-' and uppercase
    df['ID_clean'] = df[id_col].astype(str).str.split('-').str[0].str.upper()

    # Filter to 461 proteins
    df_filtered = df[df['ID_clean'].isin(proteins_461_set)].copy()
    print(f"{method}: {len(df_filtered)} proteins match 461.csv")

    # Ensure numeric columns
    df_filtered['≥80_new'] = pd.to_numeric(df_filtered['≥80_new'], errors='coerce')
    df_filtered['≥80_old'] = pd.to_numeric(df_filtered['≥80_old'], errors='coerce')
    df_filtered = df_filtered.dropna(subset=['≥80_new','≥80_old'])
    print(f"{method}: {len(df_filtered)} after dropping NaNs")

    # Count up/down/no_change per protein
    up = (df_filtered['≥80_new'] > df_filtered['≥80_old']).sum()
    down = (df_filtered['≥80_new'] < df_filtered['≥80_old']).sum()
    no_change = (df_filtered['≥80_new'] == df_filtered['≥80_old']).sum()

    summary.append({'method': method, 'up': up, 'down': down, 'no_change': no_change})

# Create summary DataFrame
summary_df = pd.DataFrame(summary)
summary_df.set_index('method', inplace=True)

# Plot stacked bar chart
ax = summary_df[['up', 'down', 'no_change']].plot(
    kind='bar',
    stacked=True,
    color=['green', 'red', 'gray'],
    figsize=(8, 5)
)

# Add labels and title
plt.ylabel("Number of proteins")
plt.xlabel("Method")
plt.title("Change in residues ≥80 for 461 proteins by method")
plt.xticks(rotation=45)
plt.legend(title="Change")
plt.tight_layout()

# Add counts on each segment for clarity
for p in ax.patches:
    height = p.get_height()
    if height > 0:
        ax.text(
            x=p.get_x() + p.get_width() / 2,
            y=p.get_y() + height / 2,
            s=int(height),
            ha='center',
            va='center',
            color='white',
            fontsize=10
        )

# Save figure
plt.savefig("stacked_bar_461_proteins.png", dpi=300)
plt.show()


# Load 461 proteins
proteins_461 = pd.read_csv("../fumigatus/461.csv")['old_id'].astype(str).str.upper().tolist()
proteins_461_set = set(proteins_461)

# Normalize ESM IDs
if 'gene_old' in esm_df.columns:
    esm_id_col = 'gene_old'
elif 'gene_name_old' in esm_df.columns:
    esm_id_col = 'gene_name_old'
else:
    esm_id_col = esm_df.columns[0]  # fallback

esm_df['ID_clean'] = esm_df[esm_id_col].astype(str).str.split('-').str[0].str.upper()

# Proteins in 461 missing from ESM
missing_in_esm = proteins_461_set - set(esm_df['ID_clean'])
print(f"Number of proteins missing in ESM: {len(missing_in_esm)}")
print("First 20 missing proteins:", list(missing_in_esm)[:20])

# Proteins present in both
present_in_esm = proteins_461_set.intersection(set(esm_df['ID_clean']))
print(f"Number of proteins present in ESM: {len(present_in_esm)}")


# Normalize the protein ID to match the cleaned IDs
protein_id = "Afu2g00090-T-p1_seed_101_sample_0"

# Filter using the normalized ID column
row = protenix1_df[protenix1_df['gene_name_old'] == protein_id]

print(row)

if not row.empty:
    with pd.option_context('display.max_columns', None):
        print(row)
else:
    print(f"{protein_id} not found in ESM DataFrame.")
