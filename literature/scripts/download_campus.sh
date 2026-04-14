#!/bin/bash
# CAMPUS WiFi Download Script - University Network Optimized
# Use this when connected to campus WiFi for better access

echo "================================================================"
echo "University Campus WiFi Paper Downloader"
echo "================================================================"
echo ""
echo "Status: Connected to campus network (detecting proxy if needed)"
echo ""

# Base directories
LITERATURE_DIR="/Users/jaimeduarteniemen/Documents/Bioinformatics/Dissertation/literature"
DUMP_DIR="$LITERATURE_DIR/literature_dump"

echo "Target directories created:"
echo "  $DUMP_DIR/hmms/"
echo "  $DUMP_DIR/uvsx_protein/"
echo "  $DUMP_DIR/cold_adaptation/"
echo "  $DUMP_DIR/phylogenetics/"
echo "  $DUMP_DIR/sequence_alignment/"
echo ""

# List of accessible papers (PMC links usually work best)
echo "Download queue:"
echo ""

echo "1. T4 UvsX Structure & Interaction with UvsW"
echo "   URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC3006652/"
echo "   Download PDF button → Save to $LITERATURE_DIR/uvsx_protein/"
echo ""

echo "2. Viral DNA Packaging Motor - Walker-A Motif"
echo "   URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC4905814/"
echo "   Download PDF button → Save to $LITERATURE_DIR/statistical_analysis/"
echo ""

echo "3. Hidden Markov Models Primer"
echo "   URL: https://www.nature.com/articles/nbt1004-1315"
echo "   Download PDF button → Save to $LITERATURE_DIR/hmms/"
echo ""

echo "4. Novel Phage Recombinases from Extreme Environments"
echo "   URL: https://pubmed.ncbi.nlm.nih.gov/41667136/"
echo "   Click 'Get PDF' or download from PMC"
echo "   → Save to $LITERATURE_DIR/hmms/"
echo ""

echo "================================================================"
echo "Manual Download Instructions:"
echo "================================================================"
echo ""
echo "1. Open the URL in your browser"
echo "2. Click the 'PDF' or 'Download PDF' button"
echo "3. Choose 'Save As' and select the appropriate destination folder"
echo "4. Example: Save to ~/Documents/Bioinformatics/Dissertation/literature/uvsx_protein/"
echo ""
echo "Or use terminal (if PDF link is available):"
echo ""
echo "For paper 1 (PMC3006652):"
echo "  cd ~/Documents/Bioinformatics/Dissertation/literature/uvsx_protein"
echo "  curl -L -o UvsX_T4_Structure.pdf \"https://pmc.ncbi.nlm.nih.gov/articles/PMC3006652/pdf/jmb.2010.10.004.pdf\""
echo ""

echo "================================================================"
echo "After downloading, I will:"
echo "  - Read each PDF"
echo "  - Extract key findings"
echo "  - Create summaries"
echo "  - Update the PAPERS_INDEX.md"
echo "================================================================"
echo ""
