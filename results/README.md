
git add .
git commit -m "Add data and results README files"
git push origin main
cd ~/telomere-fibrosis-gwas

cat > data/README.md << 'EOF'
# Data Sources

Data is not stored in this repository due to size constraints.

## Telomere Length GWAS
- Citation: Burren et al. (2024) Nature Genetics
- Accession: GCST90435144
- Available at: https://www.ebi.ac.uk/gwas/studies/GCST90435144

## Fibrosis GWAS (12 organ systems)
- Citation: Joof et al. (2026) medRxiv
- Data taken from paper

## Reference Files
- LD scores: Downloaded from Zenodo (European LD scores for LDSC)
- HapMap3 SNP list: Downloaded from Zenodo
