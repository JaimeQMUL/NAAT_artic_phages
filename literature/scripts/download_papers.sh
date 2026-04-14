#!/bin/bash
# Download papers script for literature organization
# Usage: ./download_papers.sh
# Downloads accessible PDFs to appropriate literature folders

PAPER_DIR="/Users/jaimeduarteniemen/Documents/Bioinformatics/Dissertation/literature/literature_dump"
BASE_OUTPUT="/Users/jaimeduarteniemen/Documents/Bioinformatics/Dissertation/literature"

echo "=============== Paper Downloader ============="
echo "Starting download of accessible papers..."
echo "Saving to appropriate topic folders"
echo "==============================================="
echo ""

# Function to determine category from URL
categorize_paper() {
    local url="$1"
    local title="$2"
    
    case "$url" in
        *pmc*|*pubmed*|*ovaf081*|*gkag069*)
            echo "$PAPER_DIR/hmms/$title.pdf"
            ;;
        *UvsX*|*UvsXt*|*UvsXp*|*UvsW*|*recombinase*)
            echo "$PAPER_DIR/uvsx_protein/$title.pdf"
            ;;
        *HMM*|*Markov*|*model*)
            echo "$PAPER_DIR/hmms/$title.pdf"
            ;;
        *thermal*|*temperature*|*cold*|*psychro*)
            echo "$PAPER_DIR/cold_adaptation/$title.pdf"
            ;;
        *phylogen*|*tree*|*evolution*)
            echo "$PAPER_DIR/phylogenetics/$title.pdf"
            ;;
        *alignment*|*MSA*|*sequence*|*fasta*)
            echo "$PAPER_DIR/sequence_alignment/$title.pdf"
            ;;
        *statistics*|*entropy*|*distribution*)
            echo "$PAPER_DIR/statistical_analysis/$title.pdf"
            ;;
        *visualization*|*plot*|*figure*)
            echo "$PAPER_DIR/visualization/$title.pdf"
            ;;
        *)
            echo "$PAPER_DIR/$title.pdf"
            ;;
    esac
}

# List of accessible papers to download
PAPERS=(
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC3006652/pdf/jmb.2010.10.004.pdf"
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC4905814/pdf/jmb.2016.04.029.pdf"
    "https://doi.org/10.1093/lambio/ovaf081"
    "https://www.nature.com/articles/nbt1004-1315"
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC10144138/pdf/molecules-28-03363.pdf"
)

echo "Checking university proxy/VPN status..."
echo "If you're on campus WiFi, this will have better access"
echo ""

for paper_url in "${PAPERS[@]}"; do
    echo "Processing: $paper_url"
    
    # Check if URL looks like a DOI
    if [[ "$paper_url" == *"doi.org"* ]]; then
        # Extract DOI and try to get PDF
        DOI=$(echo "$paper_url" | sed 's/.*\/doi\.org\///')
        echo "  DOI: $DOI"
        
        # Try to download PDF directly (often works from university networks)
        OUTPUT_FILE=$(categorize_paper "DOI:$DOI" "$DOI")
        
        echo "  Downloading to: $OUTPUT_FILE"
        
        # Create directory if needed
        mkdir -p "$(dirname "$OUTPUT_FILE")"
        
        # Download PDF
        if curl -L -o "$OUTPUT_FILE" "$paper_url" 2>/dev/null; then
            echo "  ✅ Downloaded successfully!"
            SIZE=$(stat -f%z "$OUTPUT_FILE" 2>/dev/null || stat -c%s "$OUTPUT_FILE" 2>/dev/null)
            echo "  Size: $SIZE bytes"
        else
            echo "  ⚠️  Download failed - may need campus WiFi or different access method"
            echo "  Try: wget $paper_url -O $OUTPUT_FILE"
        fi
    else
        echo "  Skipping - requires manual download or different approach"
    fi
    echo ""
done

echo "==============================================="
echo "Download complete!"
echo "Check folders:"
echo "  - hmms/"
echo "  - uvsx_protein/"
echo "  - cold_adaptation/"
echo "  - phylogenetics/"
echo "==============================================="
echo ""
echo "If download failed, try running from terminal:"
echo "  cd /Users/jaimeduarteniemen/Documents/Bioinformatics/Dissertation/literature"
echo "  ./scripts/download_papers.sh"
echo ""
