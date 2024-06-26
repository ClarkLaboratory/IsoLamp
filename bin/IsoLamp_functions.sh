# main script functions
function filter_reference_files() {
	
	echo "Preprocessing reads: $READS"
	
	if [ ! -d "$OUTPUT_NAME/temp_files/refs" ] 
	then
  		mkdir -p "$OUTPUT_NAME/temp_files/refs"
	fi
	# filter GTF by ensembl id
	grep $ENSG_ID $ANNOTATION > $OUTPUT_NAME/temp_files/refs/filt_chr.gtf
	ANNOTATION_FILT=$OUTPUT_NAME/temp_files/refs/filt_chr.gtf

	# get chr from 2nd line in filtered GTF to filter genome
	chr=$( cat $ANNOTATION_FILT | sed -n '2p' | awk '{ print $1}' )
	
	# check file exists and is not empty
	if [ -e "$ANNOTATION_FILT" ]
	then 	
		if [ -s "$ANNOTATION_FILT" ]
		then
			:
		else
			echo "Failed to filter GTF, check ENSG_ID is in GTF and logs file"
			exit 1
		fi
	else
		echo "Failed to filter GTF, check ENSG_ID is in GTF and logs file"
		exit 1
	fi

	# filter and index genome FASTA
	samtools faidx $GENOME $chr -o $OUTPUT_NAME/temp_files/refs/filt_chr.fa
	GENOME_FILT=$OUTPUT_NAME/temp_files/refs/filt_chr.fa

	if [ -e "$GENOME_FILT" ]
	then 	
		if [ -s "$GENOME_FILT" ]
		then
			:
		else
			echo "Failed to filter genome by chromosome, check GTF and genome have the same chromosome naming and logs file"
			exit 1
		fi
	else
		echo "Failed to filter genome by chromosome, check GTF and genome have the same chromosome naming and logs file"
		exit 1
	fi

}

function check_for_fastq_or_fasta_files() {
	
	if [ -d "$READS" ] 
	then
        	# use find to search for FASTQ or FASTA files
	        if find "$READS" -maxdepth 2 -type f \( -name "*.fastq"  -o -name "*.fastq.gz" -o -name "*.fasta" -o -name "*.fa" \) -print -quit | grep -q . 
		then
			# Check for empty files
	        	empty_files=$(find "$READS" -type f -empty)
	        	if [ -n "$empty_files" ]
			then
	           		echo "WARNING: The following empty files were found within subdirectories of $READS:"
	            		echo "$empty_files"
				echo "This may cause the program to fail."
	        	fi
		else
	        	echo "No FASTQ or FASTA files found within subdirectories of $READS"
	        	exit 1
        	fi

	else
        	echo "Directory not found: $READS"
		exit 1
	fi
    
}

function concat_files() {
	
	mkdir $OUTPUT_NAME/temp_files/reads

	# run on all subdirectories in $READS
	find "$READS" -mindepth 1 -maxdepth 1 -type d | while read -r dir
	do

		# use basename to extract just the directory name
		dir_name=$(basename "$dir")

		# gunzip zipped files
		for gzipped_file in "$READS/$dir_name"/*.gz 
		do
			if [ -f "$gzipped_file" ] 
			then
				gunzip "$gzipped_file"
			fi
		done
		
		# combine all sample read files into one file
		cat $READS/$dir_name/*.fa* > $OUTPUT_NAME/temp_files/reads/$dir_name.firstpass.fa
		
		# remove duplicated reads
		redirect_output dedupe.sh in=$OUTPUT_NAME/temp_files/reads/$dir_name.firstpass.fa out=$OUTPUT_NAME/temp_files/reads/$dir_name.fa rmn=t
		rm $OUTPUT_NAME/temp_files/reads/$dir_name.firstpass.fa

	done

	# set path to reads
	path_to_reads=$OUTPUT_NAME/temp_files/reads
    
}

function downsampling_function() {
	
	if [ -z "$downsampling" ]
	then
    		downsampling=TRUE # set to the default value if empty
  	fi

	if [ "$downsampling" == TRUE ]
	then
		
		if [ -z "$number_reads_downsample" ]
		then
    			number_reads_downsample="8000" # set to the default value if empty
  		fi

		#echo "Downsampling reads..."

		mkdir $OUTPUT_NAME/temp_files/downsampled_reads

		function downsample_for_loop() { 
			for filename in $OUTPUT_NAME/temp_files/reads/*.fa
			do

				base=$(basename "$filename")
				sample_name="${base%.*}" 
				echo "$sample_name downsampling"
				
				# downsample the reads using reformat.sh from BBMap
				reformat.sh sample=$number_reads_downsample \
				in=$OUTPUT_NAME/temp_files/reads/$sample_name.fa \
				out=$OUTPUT_NAME/temp_files/downsampled_reads/$sample_name.fa

			done
		}
		
		redirect_output downsample_for_loop

		# if statement for downsampling failing
		if [ -z "$(ls -A "$OUTPUT_NAME/temp_files/downsampled_reads")" ]
		then
   			echo "Failed to downsample reads, check logs file."
			exit 1
		else
    			:
		fi

		# set path to downsampled reads
		path_to_reads=$OUTPUT_NAME/temp_files/downsampled_reads
		
	fi

}

function mapping_genome_function() {
	
	#mkdir $OUTPUT_NAME/mapped_data

	# define max intron length for minimap2
	if [ -z "$max_intron_length" ] 
	then
		max_intron_length="400" # set to the default value if empty
	fi

	if [ -z "$splice_flank" ] 
	then
		splice_flank="yes" # set to the default value if empty
	fi
  
	echo "Mapping reads..."

    # map reads with minimap2 function
    function map_reads_minimap2() {
        local filename="$1"
        base=$(basename "$filename")
        sample_reads="${base%.*}"
        
        minimap2 -ax splice --eqx -G"${max_intron_length}"k --splice-flank="$splice_flank" $GENOME_FILT "$filename" | samtools view -bh > "$OUTPUT_NAME/mapped_data/${sample_reads}.bam" 

        # filter for primary alignments only
        samtools view -h -F 2308 "$OUTPUT_NAME/mapped_data/${sample_reads}.bam" | samtools sort - > "$OUTPUT_NAME/mapped_data/${sample_reads}_primary_sorted.bam"
    }

    # loop through each sample/barcode and map to genome with minimap2
    for filename in "$path_to_reads"/*.fa 
	do
        redirect_output map_reads_minimap2 "$filename"
    done

	# merge all BAMs 
	samtools merge -f $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam $OUTPUT_NAME/mapped_data/*sorted.bam
	samtools index $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam

	if [ -e "$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam.bai" ]
	then 
		:
	else
		echo "Failed to create BAM index, check logs file"
		exit 1
	fi
	
}

function extract_mapped_genome_reads() {

	# get list of reads mapped to genome to filter FASTQs
	samtools view $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam | awk '{print $1}' | sort | uniq > $OUTPUT_NAME/temp_files/reads_ids_mapped_to_genome.txt
	
	if [ -z "$primer_site_based_filter" ]
	then
    		primer_site_based_filter=FALSE # set to the default value if empty
  	fi

	# extract primers per line
	if [ "$primer_site_based_filter" == TRUE ]
	then

		mkdir $OUTPUT_NAME/temp_files/primers
		mkdir $OUTPUT_NAME/temp_files/primers/forward
		mkdir $OUTPUT_NAME/temp_files/primers/reverse

		# forward
		counts_lines_forward=$(awk 'END{print NR}' $forward_primers)
		for i in $(seq 1 $counts_lines_forward)
		do
			sed -n "${i}p" $forward_primers > $OUTPUT_NAME/temp_files/primers/forward/forward_${i}.bed
		done

		# reverse
		counts_lines_reverse=$(awk 'END{print NR}' $reverse_primers)
		for i in $(seq 1 $counts_lines_reverse)
		do
			sed -n "${i}p" $reverse_primers > $OUTPUT_NAME/temp_files/primers/reverse/reverse_${i}.bed
		done

		# get read start and end windows

		cat $forward_primers $reverse_primers > $OUTPUT_NAME/temp_files/primers/combined.bed
		lines_in_combined_bed=$(cat $OUTPUT_NAME/temp_files/primers/combined.bed | wc -l)

		if [ $((counts_lines_forward + counts_lines_reverse)) -eq $lines_in_combined_bed ]
		then
			:
		else
			echo "Error in primer filtering. Hit 'Enter' key after the last line in your bed files to ensure the lines are read correctly."
			exit 1
		fi

		reads_start=$(sort -k2 -n $OUTPUT_NAME/temp_files/primers/combined.bed | head -1 | awk '{print $2}')
		window_start=$(($reads_start - 100))
		
		reads_end=$(sort -k3 -n $OUTPUT_NAME/temp_files/primers/combined.bed | tail -1 | awk '{print $3}')
		window_end=$(($reads_end + 100))

		bedtools bamtobed -i $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam > $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bed

		awk -v a="$window_start" -v b="$window_end" '$2 > a && $3 < b' $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bed | awk '{print $4}' > $OUTPUT_NAME/mapped_data/reads_in_window.txt

		samtools view -bh -N $OUTPUT_NAME/mapped_data/reads_in_window.txt $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam > $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged_in_window.bam

		# forward
		forward_primer_count=$(ls -1 "$OUTPUT_NAME/temp_files/primers/forward" | wc -l)
		for i in $(seq 1 $forward_primer_count)
		do
		bedtools intersect -split -wa -u -F 0.8 -abam $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged_in_window.bam -b $OUTPUT_NAME/temp_files/primers/forward/forward_${i}.bed > $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_overlap_forward_primer_${i}.bam
		done

		samtools merge -o $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_merged_overlap_forward_primers.bam $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_overlap_forward_primer*.bam

		# reverse
		reverse_primer_count=$(ls -1 "$OUTPUT_NAME/temp_files/primers/reverse" | wc -l)
		for i in $(seq 1 $reverse_primer_count)
		do
		bedtools intersect -split -wa -u -F 0.8 -a $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_merged_overlap_forward_primers.bam -b $OUTPUT_NAME/temp_files/primers/reverse/reverse_${i}.bed > $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_overlap_reverse_primer_${i}.bam
		done

		samtools merge -o $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_merged_overlap_reverse_primers.bam $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_overlap_reverse_primer*.bam

		# get list of reads mapped to primers and genome to filter FASTQs
		samtools view $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_merged_overlap_reverse_primers.bam | awk '{print $1}' | sort | uniq > $OUTPUT_NAME/temp_files/reads_ids_full_length.txt

		# get reads in primer filter
		list_for_filtering_reads=$OUTPUT_NAME/temp_files/reads_ids_full_length.txt
		
		

	else
		# get reads mapped to genome
		list_for_filtering_reads=$OUTPUT_NAME/temp_files/reads_ids_mapped_to_genome.txt
	fi

	mkdir $OUTPUT_NAME/temp_files/mapped_reads
	
	function extract_mapped_reads_loop() {
		for filename in $path_to_reads/*.fa
		do
			base=$(basename $filename)
			sample_reads="${base%.*}"
			
			filterbyname.sh in=$path_to_reads/${sample_reads}.fa out=$OUTPUT_NAME/temp_files/mapped_reads/${sample_reads}.fa names=$list_for_filtering_reads include=t
		done
	}

	redirect_output extract_mapped_reads_loop
	path_to_reads=$OUTPUT_NAME/temp_files/mapped_reads

	# print total reads per barcode mapped to genome (and full length if primer filtering is TRUE)
	for filename in $OUTPUT_NAME/temp_files/mapped_reads/*.fa
	do
		base=$(basename "$filename")
		sample_name="${base%.*}" 
				
		cat $OUTPUT_NAME/temp_files/mapped_reads/$sample_name.fa | grep "^[>@]" | wc -l >> $OUTPUT_NAME/temp_files/reads_per_barcode_mapped_genome.txt 
	done

}

function get_accurate_reads() {
	
	if [ -z "$extract_high_accuracy_reads" ]
	then
    		extract_high_accuracy_reads=TRUE # set to the default value if empty
  	fi

	if [ "$extract_high_accuracy_reads" == TRUE ]
	then
		
		if [ -z "$minimum_read_accuracy" ]
		then
    			minimum_read_accuracy=0.95 # set to the default value if empty
  		fi

		#echo "Extracting high accuracy reads..."
		if [ "$primer_site_based_filter" == TRUE ]
		then
			redirect_output Rscript $SCRIPT_DIR/bin/read_list_high_accuracy.R -i $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_merged_overlap_reverse_primers.bam -a $minimum_read_accuracy -o $OUTPUT_NAME
			samtools view -N $OUTPUT_NAME/temp_files/reads_above_accuracy_minimum.txt -o $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_accuracy_reads.bam $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_merged_overlap_reverse_primers.bam
		else
			redirect_output Rscript $SCRIPT_DIR/bin/read_list_high_accuracy.R -i $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam -a $minimum_read_accuracy -o $OUTPUT_NAME
			samtools view -N $OUTPUT_NAME/temp_files/reads_above_accuracy_minimum.txt -o $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_accuracy_reads.bam $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam
		fi
		
		samtools index $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_accuracy_reads.bam

		bam_for_bambu=$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_accuracy_reads.bam

		# Check completed
		if [ -e "$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_accuracy_reads.bam" ]
		then
			if [ -s "$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_accuracy_reads.bam" ]
			then
				:
			else
				echo "Failed to extract high accuracy reads, check logs file"
				exit 1
			fi
		else
			echo "Failed to extract high accuracy reads, check logs file"
			exit 1
		fi

	fi

}

function get_JWRs() {
	
	if [ -z "$extract_high_quality_SJs" ]
	then
    		extract_high_quality_SJs=FALSE # set to the default value if empty
  	fi

	if [ "$extract_high_quality_SJs" == TRUE ] && [ "$extract_high_accuracy_reads" == TRUE ]
	then
		
		if [ -z "$junction_window" ]
		then
    			junction_window=20 # set to the default value if empty
  		fi

		if [ -z "$JAQ" ]
		then
    			JAQ=0.9 # set to the default value if empty
  		fi
		
		# update BAM

		python3 $SCRIPT_DIR/bin/NanoSplicer_JWR_checker.py --window $junction_window $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_accuracy_reads.bam $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_JAQ.csv
		Rscript $SCRIPT_DIR/bin/read_list_high_JAQ.R -i $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_JAQ.csv -j $JAQ -o $OUTPUT_NAME
		samtools view -N $OUTPUT_NAME/temp_files/reads_above_JAQ_minimum.txt -o $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_JAQ_reads.bam $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_accuracy_reads.bam

		

		bam_for_bambu=$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_JAQ_reads.bam

		# Check JWR completed
		if [ -e "$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_JAQ_reads.bam" ]
		then
			if [ -s "$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_JAQ_reads.bam" ]
			then
				:
			else
				echo "Failed to extract splice junctions, check logs file"
				exit 1
			fi
		else
			echo "Failed to extract splice junctions, check logs file"
			exit 1
		fi
	
	elif [ "$extract_high_quality_SJs" == TRUE ] && [ "$extract_high_accuracy_reads" == FALSE ]
	then
    	echo "Cannot extract high quality splice junctions if high accuracy reads are not extracted first. Set extract_high_accuracy_reads=TRUE."
		echo "Skipping SJ extraction."
		bam_for_bambu=$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam
	elif [ "$extract_high_quality_SJs" == FALSE ] && [ "$extract_high_accuracy_reads" == FALSE ]
	then
		bam_for_bambu=$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam
	fi

}

function set_skip_mapping_vars() {
	
	if ! test -f $OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam
	then
		echo "No BAM file found, run mapping first."
		exit 1
	fi

	# define
	ANNOTATION_FILT=$OUTPUT_NAME/temp_files/refs/filt_chr.gtf
	GENOME_FILT=$OUTPUT_NAME/temp_files/refs/filt_chr.fa

	if [ "$extract_high_accuracy_reads" == TRUE ] && [ "$extract_high_quality_SJs" == TRUE ] 
	then
		bam_for_bambu=$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_JAQ_reads.bam
	elif [ "$extract_high_accuracy_reads" == TRUE ] && [ "$extract_high_quality_SJs" == FALSE ]
	then
		bam_for_bambu=$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_filtered_high_accuracy_reads.bam
	else
		bam_for_bambu=$OUTPUT_NAME/mapped_data/${OUTPUT_NAME}_primary_merged.bam
	fi

 	if [ "$primer_site_based_filter" == TRUE ]
	then
		forward_primer_count=$(ls -1 "$OUTPUT_NAME/temp_files/primers/forward" | wc -l)
		reverse_primer_count=$(ls -1 "$OUTPUT_NAME/temp_files/primers/reverse" | wc -l)
  	fi
   
	num_of_barcodes=$( ls -lh $OUTPUT_NAME/temp_files/reads/*.f* | wc -l )
	path_to_reads=$OUTPUT_NAME/temp_files/mapped_reads
	
}

function run_bambu_function() {
	
	echo "Identifying isoforms with Bambu..."

	if [ -z "$bambu_ndr" ] 
  	then
    		bambu_ndr=1 # set to the default value
  	fi

	if [ -z "$bambu_min_gene_fraction" ] 
  	then
    		bambu_min_gene_fraction=0.001 # set to the default value
  	fi

	if [ ! -d "$OUTPUT_NAME/bambu" ] 
	then
  		mkdir -p "$OUTPUT_NAME/bambu"
	fi


	# run bambu Rscript 
	redirect_output Rscript $SCRIPT_DIR/bin/bambu_tx_discovery.R -b $bam_for_bambu \
		-f $GENOME_FILT \
		-t $ANNOTATION_FILT \
		-n $bambu_ndr \
		-g $bambu_min_gene_fraction \
		-o $OUTPUT_NAME/bambu

	if [ -e "$OUTPUT_NAME/bambu/extended_annotations.gtf" ]
	then
		if [ -s "$OUTPUT_NAME/bambu/extended_annotations.gtf" ]
		then
			:
		else
			echo "Failed to create Bambu annotations, check logs file"
			exit 1
		fi
	else
		echo "Failed to create Bambu annotations, check logs file"
		exit 1
	fi

}

function read_count_function_first_pass() {
	
	if ls $OUTPUT_NAME/temp_files/${OUTPUT_NAME}*.txt 1> /dev/null 2>&1 
	then
		rm $OUTPUT_NAME/temp_files/${OUTPUT_NAME}*.txt
	fi

	if ls $OUTPUT_NAME/temp_files/*.gtf 1> /dev/null 2>&1 
	then
		rm $OUTPUT_NAME/temp_files/*.gtf
	fi

	# filter for bambu isoforms with read count > 2 to remove very lowly expressed isoforms
	cat "$OUTPUT_NAME/bambu/counts_transcript.txt" | awk '{ if ($3 > 2 ) print $1 }' | tail -n +2 > "$OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_read_count_min_first_pass.txt"
	
	# subset bambu GTF for the isoforms that passed threshold
	cat $OUTPUT_NAME/bambu/extended_annotations.gtf | grep -wf $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_read_count_min_first_pass.txt > $OUTPUT_NAME/bambu/extended_annotations_first_pass.gtf
	
	# replace some Bambu naming conventions
	cat $OUTPUT_NAME/bambu/extended_annotations_first_pass.gtf | sed 's/tx./tx/g' | sed 's/BambuTx/tx/g' > $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_first_pass.gtf
	cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_read_count_min_first_pass.txt | sed 's/tx./tx/g' | sed 's/BambuTx/tx/g' | awk '{ print $1"\t"$3}' > $OUTPUT_NAME/temp_files/updated_transcriptome_first_pass.txt

	# check file exists and is not empty
	if [ -e "$OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_first_pass.gtf" ]
	then 	
		if [ -s "$OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_first_pass.gtf" ]
		then
			:
		else
			echo "Read count threshold first pass failed, check logs file"
			exit 1
		fi
	else
		echo "Read count threshold first pass failed, check logs file"
		exit 1
	fi

}

function primer_site_function() {

	if [ "$primer_site_based_filter" == TRUE ]
	then

		#echo "Filtering based on primer BED files..."
		
		# intersect forward and reverse primers with GTF
		# forward
		for i in $(seq 1 $forward_primer_count)
		do
		cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_first_pass.gtf | awk '{ if ($3 == "exon") print }' | bedtools intersect -wa -u -F 0.8 -a - -b $OUTPUT_NAME/temp_files/primers/forward/forward_${i}.bed | grep -o 'transcript_id "[^"]*' | cut -d' ' -f2 | sed 's/"//g' > $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_forward_primer_exons_${i}.txt
		done

		cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_forward_primer_exons_*.txt > $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_forward_primer_exons_all.txt

		# reverse
		for j in $(seq 1 $reverse_primer_count)
		do
		cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_first_pass.gtf | awk '{ if ($3 == "exon") print }' | bedtools intersect -wa -u -F 0.8 -a - -b $OUTPUT_NAME/temp_files/primers/reverse/reverse_${j}.bed | grep -o 'transcript_id "[^"]*' | cut -d' ' -f2 | sed 's/"//g' > ${OUTPUT_NAME}/temp_files/${OUTPUT_NAME}_isoforms_reverse_primer_exons_${j}.txt
		done

		cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_reverse_primer_exons_*.txt > $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_reverse_primer_exons_all.txt
		
		# combine isoforms found in both forward and reverse primers
		comm -12 <(sort $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_forward_primer_exons_all.txt) <(sort $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_reverse_primer_exons_all.txt) > $OUTPUT_NAME/temp_files/updated_transcriptome.txt
		
		# subset GTF that passed counts filtering with new list of isoforms that match primers
		cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_first_pass.gtf | grep -wf $OUTPUT_NAME/temp_files/updated_transcriptome.txt > $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_primer.gtf

		if [ -e "$OUTPUT_NAME/temp_files/updated_transcriptome.txt" ]
		then 	
			if [ -s "$OUTPUT_NAME/temp_files/updated_transcriptome.txt" ]
			then
				:
			else
				echo "Failed to filter based on primers, check logs file"
				exit 1
			fi
		else
			echo "Failed to filter based on primers, check logs file"
			exit 1
		fi

		# sum transcript counts remaining after primer filtering
		total_kept_counts=$(cat "$OUTPUT_NAME/bambu/counts_transcript.txt" | grep -wf "$OUTPUT_NAME/temp_files/updated_transcriptome.txt" | awk '{total=total+$3} END{printf "%.0f", total}')
		# sum transcript counts from Bambu
		total_bambu_counts=$(cat "$OUTPUT_NAME/bambu/counts_transcript.txt" | awk '{total=total+$3} END{printf "%.0f", total}')
		
		# calculate 90% of the total transcript counts
		primer_threshold=$(printf "%.0f" "$(echo "$total_bambu_counts * 0.7" | bc)")

		# check that at least 90% of the counts are remaining after primer filtering
		if [ "$((total_kept_counts))" -lt "$((primer_threshold))" ] 
		then
			echo "WARNING: primer filtering removed isoforms with >30% of counts assigned by Bambu"
		fi

	elif [ "$primer_site_based_filter" == FALSE ]
	then
		# use isoforms that passed read count filtering 
		cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_first_pass.gtf > $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_primer.gtf
		cat $OUTPUT_NAME/temp_files/updated_transcriptome_first_pass.txt > $OUTPUT_NAME/temp_files/updated_transcriptome.txt
	fi

}

function create_metatranscriptome() {
	
	#echo "Mapping to updated transcriptome with novel isoforms..."
	if [ ! -d "$OUTPUT_NAME/updated_transcriptome" ] 
	then
  		mkdir -p "$OUTPUT_NAME/updated_transcriptome"
	fi

	# create an updated transcriptome fasta with gffread using the bambu isoforms
	redirect_output gffread -w $OUTPUT_NAME/updated_transcriptome/updated_transcriptome.fa -g $GENOME_FILT $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_primer.gtf
	
	if [ -e "$OUTPUT_NAME/updated_transcriptome/updated_transcriptome.fa" ]
  	then 
		:
 	else
		echo "Failed to create updated transcriptome with gffread, check logs file"
		exit 1
  	fi

	# create a salmon index of the updated transcriptome
	redirect_output salmon index -t $OUTPUT_NAME/updated_transcriptome/updated_transcriptome.fa -i $OUTPUT_NAME/updated_transcriptome -k 31
	
	if [ -e "$OUTPUT_NAME/updated_transcriptome/versionInfo.json" ]
  	then 
		:
 	else
		echo "Failed to index with salmon, check logs file"
		exit 1
  	fi

}

function remapping_function() {
	
	if [ ! -d "$OUTPUT_NAME/updated_transcriptome/salmon_quants" ] 
	then
  		mkdir -p "$OUTPUT_NAME/updated_transcriptome/salmon_quants"
	fi

	# function to loop through reads from each sample/barcode and quantify with salmon
	running_salmon() {
		for filename in $path_to_reads/*.fa
		do

			base=$(basename $filename)
			sample_reads="${base%.*}"
			
			salmon quant -i $OUTPUT_NAME/updated_transcriptome \
			-l A \
			-r $path_to_reads/${sample_reads}.fa \
			-o $OUTPUT_NAME/updated_transcriptome/salmon_quants/${sample_reads} \
			-z $OUTPUT_NAME/updated_transcriptome/salmon_quants/${sample_reads}/mapped_${sample_reads}.sam \
			--auxDir aux \
			--writeUnmappedNames
		
		done
	}

	redirect_output running_salmon
		
	if [ -d "$OUTPUT_NAME/updated_transcriptome/salmon_quants" ]
  	then
		if [ "$(ls -A $OUTPUT_NAME/updated_transcriptome/salmon_quants)" ] 
		then
     			:
		else
    			echo "Failed to quantify with salmon, check logs file"
			exit 1
		fi
 	else
		echo "Failed to quantify with salmon, check logs file"
		exit 1
  	fi

}

function combine_quants_function() {
	
	if [ -z "$TPM_minimum" ]
	then
    		TPM_minimum="5000" # set to the default value if empty
  	fi

	if [ -z "$samples_minimum" ]
	then
		samples_minimum=$(echo "scale=0; $num_of_barcodes/4" | bc -l) # set to the default value if empty
  	fi

	if (( $samples_minimum == 0 ))
	then
		samples_minimum=1
  	fi
	
	echo "Generating output files..."

	# run Rscript to combine quant.sf files produced by salmon to produce isoform counts
	Rscript $SCRIPT_DIR/bin/combine_salmon_quants.R -e $ENSG_ID -m $TPM_minimum -s $samples_minimum -o $OUTPUT_NAME

	

}

function read_count_function_second_pass() {

	# create list in text file of isoforms that passed threshold
	# command uses _ and , as sep, prints first column of tx names, removes the word Bambu in tx names, and removes the column header
	cat $OUTPUT_NAME/${OUTPUT_NAME}_counts.csv | awk -F '[,]' '{print $2}' | sed 's/Bambu//g' | tail -n +2 > $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_read_count_min_second_pass.txt
	
	# filter GTF based on isoforms 
	if ls ${OUTPUT_NAME}/*.gtf 1> /dev/null 2>&1 
	then
		rm ${OUTPUT_NAME}/*.gtf
	fi
	cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_primer.gtf | grep -wif $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_read_count_min_second_pass.txt > $OUTPUT_NAME/${OUTPUT_NAME}_isoforms.gtf

	
	if [ -e "$OUTPUT_NAME/${OUTPUT_NAME}_isoforms.gtf" ]
	then 	
		if [ -s "$OUTPUT_NAME/${OUTPUT_NAME}_isoforms.gtf" ]
		then
			:
		else
			echo "Read count threshold second pass failed, check logs file"
			exit 1
		fi
	else
		echo "Read count threshold second pass failed, check logs file"
		exit 1
	fi

}

function gffcomp_function() {
	if [ ! -d "$OUTPUT_NAME/annotated_isoforms" ] 
	then
  		mkdir -p "$OUTPUT_NAME/annotated_isoforms"
	fi

	# create file of bambu tx classes for remaining isoforms
	{ head -n 1 $OUTPUT_NAME/bambu/bambu_tx_classes.txt && grep -wif $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_read_count_min_second_pass.txt $OUTPUT_NAME/bambu/bambu_tx_classes.txt; } > $OUTPUT_NAME/annotated_isoforms/${OUTPUT_NAME}_bambu_isoform_classes.txt
	
	# run gffcompare on the new isoforms and original reference GTF
	redirect_output gffcompare -r $ANNOTATION_FILT -o $OUTPUT_NAME/annotated_isoforms/gffcomp $OUTPUT_NAME/${OUTPUT_NAME}_isoforms.gtf
	
	# put gffcompare outputs in directory
	mv $OUTPUT_NAME/gffcomp.* $OUTPUT_NAME/annotated_isoforms

}

function plot_PCA() {

	if [ -z "$grouping_data" ] 
	then
		grouping_data=NULL # set to the default value if empty
	fi

	redirect_output Rscript $SCRIPT_DIR/bin/plot_PCA.R -i $OUTPUT_NAME/$OUTPUT_NAME"_counts.csv" -o "$OUTPUT_NAME/$OUTPUT_NAME" -g $grouping_data

}

function plot_accuracy () {
	
	redirect_output Rscript $SCRIPT_DIR/bin/accuracy.R $OUTPUT_NAME/mapped_data/$OUTPUT_NAME"_primary_merged.bam" "$OUTPUT_NAME/$OUTPUT_NAME"

}

function compare_proportions_test() {
	
	# only run if grouping is provided
	if [ -z "$grouping_data" ] 
	then
		:
	else
		#echo "Comparing isoform proportions between groups..."
		redirect_output Rscript $SCRIPT_DIR/bin/prop_t_test.R -i $OUTPUT_NAME/$OUTPUT_NAME"_proportions.csv" -o $OUTPUT_NAME -g $grouping_data
	fi

}

function generate_report_function() {

	if [ -z "$extract_high_accuracy_reads" ]
	then
    		extract_high_accuracy_reads=TRUE # set to the default value if empty
			minimum_read_accuracy=0.95
  	fi

	if [ "$extract_high_accuracy_reads" == TRUE ]
	then
		number_reads_acc=$(cat $OUTPUT_NAME/temp_files/reads_above_accuracy_minimum.txt | sort | uniq | wc -l)
	fi

	if [ -z "$extract_high_quality_SJs" ]
	then
    		extract_high_quality_SJs=FALSE # set to the default value if empty
  	fi

	if [ "$extract_high_quality_SJs" == TRUE ]
	then
		number_reads_SJ=$(cat $OUTPUT_NAME/temp_files/reads_above_JAQ_minimum.txt | sort | uniq | wc -l)
	fi

	# report calculations
	num_of_barcodes=$( ls -lh $OUTPUT_NAME/temp_files/reads/*.f* | wc -l )
	
	if ls $OUTPUT_NAME/temp_files/reads/*.gz 1> /dev/null 2>&1 
	then
		:
	else
		total_reads_input=$( cat $OUTPUT_NAME/temp_files/reads/*.f* | grep "^[>@]" | wc -l )
	fi

	if ls $OUTPUT_NAME/temp_files/downsampled_reads/*.gz 1> /dev/null 2>&1 
	then
		:
	else
		total_reads_post_downsample=$(cat $OUTPUT_NAME/temp_files/downsampled_reads/*.f* | grep "^[>@]" | sort | uniq | wc -l)
	fi
	
	number_reads_mapped=$(cat $OUTPUT_NAME/temp_files/reads_ids_mapped_to_genome.txt | wc -l)
	number_reads_inside_window=$(cat $OUTPUT_NAME/mapped_data/reads_in_window.txt | wc -l)
	number_reads_full_length=$(cat $OUTPUT_NAME/temp_files/reads_ids_full_length.txt | wc -l)
	filtered_transcripts_known=$( cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_read_count_min_second_pass.txt | grep -vi tx | wc -l )
	filtered_transcripts_novel=$( cat $OUTPUT_NAME/temp_files/${OUTPUT_NAME}_isoforms_read_count_min_second_pass.txt | grep -i tx | wc -l )
	reads_remaining_final=$(cat $OUTPUT_NAME/temp_files/remaining_read_sum.txt)

cat << EOT >> $OUTPUT_NAME/${OUTPUT_NAME}_report.txt
`date "+%Y-%m-%d %H:%M:%S"`

Filters applied:
	Downsampling: $downsampling
	Downsampled to number of reads: $number_reads_downsample
	Primer site filter: $primer_site_based_filter
	Extract high accuracy reads: $extract_high_accuracy_reads
	Accuracy minimum: $minimum_read_accuracy
	Extract high quality splice junctions: $extract_high_quality_SJs
	Splice junction accuracy minimum: $JAQ
	TPM minimum: $TPM_minimum
	Samples/barcodes required to meet TPM minimum: $samples_minimum
	Bambu novel discovery rate (NDR): $bambu_ndr
	Bambu minimum read fraction by gene: $bambu_min_gene_fraction

Number of input samples/barcodes: $num_of_barcodes
Total number of reads in samples/barcodes: $total_reads_input

Total number of reads in samples/barcodes post-downsampling: $total_reads_post_downsample
Number of reads mapped to genome: $number_reads_mapped
Number of reads mapped inside primer window: $number_reads_inside_window
Number of reads mapped to genome and overlap forward and reverse primers: $number_reads_full_length
Number of high accuracy reads used to identify isoforms: $number_reads_acc
Number of reads with accurate SJs used to identify isoforms: $number_reads_SJ

Known isoforms identified: $filtered_transcripts_known
Novel isoforms identified: $filtered_transcripts_novel

Number of reads mapped to known/novel isoforms: $reads_remaining_final

EOT
}
