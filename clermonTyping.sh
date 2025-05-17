#!/usr/bin/env bash

#################################################
########## Clermont Typing Pipeline #############
#################################################
# From a set of contigs in FASTA format:
# 1] Launch Mash for identifying phylogroup
# 2] Make a BLAST database
# 3] Launch BLAST on a primers FASTA file
# 4] Launch in silico PCR for obtaining phylogroup
# 5] Generating reports
#
# Current "OLD" version: 24.02 (February 2024)
# Dantas Lab Release: 1.0.0 (May 2025)
version="Clermont Typing v1.0.0 (May 2025)"

# Original Contact: antoine.bridier-nahmias@inserm.fr
# Current Contact for Dantas Lab Release: caelanjmiller@wustl.edu

# More dynamic method of PATH searching for identifying "hard-coded" bin & data directories
SPACK_INSTALLATION_DIR=$(dirname $(dirname $(realpath $0)))
DATA_DIR="${SPACK_INSTALLATION_DIR}/data"
BIN_DIR="${SPACK_INSTALLATION_DIR}/bin"

# Default threshold = 0 (disabled)
THRESHOLD=0
# Default name = date
DATE=$( date "+%F_%H%M%S")
NAME=analysis_$DATE
# BLAST settings
PRIMERS="${DATA_DIR}/primers.fasta"
PERC_IDENTITY=90
BLAST_TASK='blastn'
# Mash settings
DB_MASH="${DATA_DIR}/mash/mash_reference.msh"
# Flag for minimal
MINIMAL=0

# Check if primers.fasta exists
if [ ! -f "${PRIMERS}" ]; then
  echo "Error: primers.fasta not found"
  exit 1
fi

# Check if the Mash database exists
if [ ! -f "${DB_MASH}" ]; then
  echo "Error: mash_reference.msh not found"
  exit 1
fi


function usage(){
	printf "Script usage :\n"
	printf "\t-h					: Print this message and exit\n"
	printf "\t-v					: Print the version and exit\n"
	printf "\t--fasta					: FASTA contigs file(s). If multiple files, they must be separated by an arobase (@) value\n"
	printf "\t--name					: Name for this analysis (optional)\n"
	printf "\t--threshold				: Option for Clermont Typing, do not use contigs under this size (optional)\n"
	printf "\t--minimal				: Output a minimal set of files (optional)\n"
	printf "\t--fastafile				: File with path of FASTA contig file.  One file by line (optional)\n"
	printf "\t--summary				: File with path of *_phylogroups.txt. One file by line (optional)\n"
}

function mash_analysis(){
	echo "===== Running Mash ====="
	"${BIN_DIR}/mash" screen -w "$DB_MASH" "$FASTA" > "$WORKING_DIR/${FASTA_NAME}_mash_screen.tab"
}

function blast_analysis(){
	echo "===== Creating BLAST database ====="
	echo "makeblastdb -in $FASTA -input_type fasta -out $WORKING_DIR/db/$NAME -dbtype nucl"
	makeblastdb -in "$FASTA" -input_type fasta -out "$WORKING_DIR/db/$FASTA_NAME" -dbtype nucl
	if [ $? -eq 0 ]
	then
        echo "===== Running BLAST ====="
        blastn -query "$PRIMERS" -perc_identity "$PERC_IDENTITY" -task "$BLAST_TASK" -word_size 6 -outfmt 5 -db "$WORKING_DIR/db/$FASTA_NAME" -out "$WORKING_DIR/$FASTA_NAME.xml"
        error=0
	else
        echo "Error Detected! Exiting pipeline..."
        error=1
	fi

}

function report_calling(){
	# Rscript = path to clermontReport.R
	# clermont_out = path to clermonTyping output 
	# namus = report name
	# out_dir = self explanatory!
	echo "============= Generating Report ==============="
	rscript=$1
	shift
	clermont_out=$1
	shift
	namus=$1
	shift
	out_dir=$1

	modif_script="${out_dir}/${namus}.R"
	cp "${rscript}" "${modif_script}"

	sed -i "s:TARTAMPION:$clermont_out:g" "${modif_script}"

	Rscript --slave -e "library(markdown); sink('/dev/null');rmarkdown::render('${modif_script}')"
}

function add_mash_group() {
    in_file=$1
    Rscript "${BIN_DIR}/add_mash_minimal.R" ${in_file}
}

if [ $# == 0 ]
then
	 usage
	 exit 1
fi

while [[ $# -gt 0 ]]
do
    case "$1" in
    -v)
        echo "$version"
        usage
        exit 0
        ;;
    -h)
        usage
        exit 0
        ;;
    --fastafile) 
        FASTA_FILE="$2";
        shift
        ;;    
    --fasta) 
        FASTAS="$2";
        shift
        ;;
    --name) 
        NAME="$2";
        shift
        ;;
    --minimal)
        MINIMAL=1;
        ;;	
    --threshold) 
        THRESHOLD="$2";
        shift
        ;;
    --summary) 
        SUMMARY="$2";
        shift
        ;;    

    --) shift; break;;
    esac
    shift
done


if [ -z "$FASTAS" ] && [ -z "$FASTA_FILE" ] && [ -z "$SUMMARY" ] 
then
    echo "Missing the contigs file. Option --fasta or --fastafile"
    usage
    exit 1
elif [ -n "$FASTAS" ] && [ -n "$FASTA_FILE" ]
then
    echo "Too many parameters. Option --fasta or --fastafile"
    usage  
    exit 1
elif [ -n "$SUMMARY" ] 
then
    NAME='Summary'
    echo "You asked for a Clermont Typing analysis named $NAME of phylogroups."
    if  [ -n "$FASTA_FILE" ] || [ -n "$FASTAS" ] 
    then
        echo "Too many parameters. Option --fasta or --fastafile, or --summary"
        usage
        exit 1
    fi

else
    echo "You asked for a Clermont Typing analysis named $NAME of phylogroups on $FASTAS with a minimum contig size of $THRESHOLD."
fi

if [ ! -d $NAME ]
then
	mkdir $NAME
fi
CURRENT_DIR=$(pwd)
WORKING_DIR=$CURRENT_DIR/$NAME


if [ -z "$SUMMARY" ]
then
    declare -a LIST_FILES
    if [ -z "$FASTA_FILE" ]
    then
        IFS='@' read -ra ARRAY_FASTA <<< "$FASTAS"
        LIST_FILES="${ARRAY_FASTA[@]}"
    else    
        while IFS=$'\n' read -r line 
        do 
            sample=${line%\\n}  
            LIST_FILES+=("${sample}")
        done < "${FASTA_FILE}"
    fi

    # Analysis of each FASTA file
    for FASTA in ${LIST_FILES[@]}
    do
        if [ -f "$FASTA" ] && [ -n "$FASTA" ]
        then
            # Rename File
            BASE_NAME_FASTA=$(basename "$FASTA")
            FASTA_NAME=${BASE_NAME_FASTA%£*}
            echo "============== Analysis of ${FASTA_NAME} =================="
            cp "$FASTA" "$WORKING_DIR/${FASTA_NAME}"
            FASTA="$WORKING_DIR/$FASTA_NAME"
            
            ##### Step 1: Mash Analysis #####
            # Generate ${FASTA_NAME}_mash_screen.tab
            mash_analysis
            
            ##### Step 2: BLAST #############
            # Generate ${FASTA_NAME}.xml
            blast_analysis
            if [ $error -gt 0 ]
            then
                printf "%s\t\t\t\tNA\t%s__mash_screen.tab\n" "$FASTA_NAME" "$FASTA_NAME" >> "$WORKING_DIR/${NAME}_phylogroups.txt"
            else
                ##### Step 3: Clermon Typing #####
                echo "====== Clermon Typing ====="
                results=$("${BIN_DIR}/clermont.py" -x "${WORKING_DIR}/${FASTA_NAME}.xml" -s "$THRESHOLD")
                printf "%s\t%s\t%s_mash_screen.tab\n" "$FASTA_NAME" "$results" "$FASTA_NAME" >> "$WORKING_DIR/${NAME}_phylogroups.txt"
            fi
        else
            echo "$FASTA doesn't exist"
        fi
    done
else
    while IFS=$'\n' read -r line 
    do 
        IFS='/' read -ra ARRAY_PATH <<< "$line"
        DIR="${ARRAY_PATH[0]}"
    
        IFS=$'\n'
        for line2 in $(cat $line)
        do
            IFS=$'\t'
            i=0
            for line3 in $line2
            do
                if (( i < 5 ))
                then
                    printf "%s\t" "$line3" >> "${WORKING_DIR}/${NAME}_phylogroups.txt"
                else
                    printf "../%s/%s" "$DIR" "$line3" >> "${WORKING_DIR}/${NAME}_phylogroups.txt"
                fi
                ((i++))
            done
        done
        printf "\n" >> "${WORKING_DIR}/${NAME}_phylogroups.txt"
    done < "${SUMMARY}"
fi


##### Step 4: Reporting #########
if  [ ${MINIMAL} -eq 0 ] 
then
    report_calling "${BIN_DIR}/clermontReport.R" "$WORKING_DIR/${NAME}_phylogroups.txt" "$NAME" "$WORKING_DIR" add_mash_group
    add_mash_group "$WORKING_DIR/${NAME}_phylogroups.txt 0.95"
else 
    add_mash_group "$WORKING_DIR/${NAME}_phylogroups.txt 0.95"
fi

echo "============== End =================="
exit 0




