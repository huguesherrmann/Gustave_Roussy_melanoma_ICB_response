#!/bin/bash

########## -1)  input variables ################
###########################################

while getopts "b:k:o:l:t:f:" opt
do
        case $opt in
                b)      metabam_file=$OPTARG
                        ;;
                k)      kmer_file=$OPTARG
                        ;;
                o)      output_dir=$OPTARG
                        ;;
                l)      lib_type=$OPTARG
                        ;;
                t)      samtools_threads=$OPTARG
                        ;;
                f)      fract=$OPTARG
                        ;;
                *)      echo -e ${RED}"ERROR: invalid option provided"${NOCOLOR}
                        ;;
        esac
done


#column of kmers (tag) in the table
col_kmer=1

#column of contigs in the table
col_contig=2

#awk script to determine the strand from direct and indirect grep on bam files (strand + & -)
#indirect grep is a grep after the reverse complementation of the sequence (here the sequence of the kmer/tag)
#if empty, we assume getFixedKmerStrand.awk is in the same directory as getStrandedKmers.sh
getFixedKmerStrand="/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_k_mer_expression/scripts/get_fixed_kmers_strand.awk"

#kmer_file="/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/characterize_stable_regulons/gr1234/annotation_and_sequences/nr_stable_intergenic_sequences.tsv"

#output directory (if it doesn't exist yet, it will be automatically created)
#output_dir="/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/bi_direction_transcription/gr1234/D437T40/"

#if you already have a metabam file, put the full path to the file here
#metabam_file="/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/star/gr1234/bam/D437T40/D437T40.bam"

#samtools location
samtools="samtools"

#chose the library type ; you can check your bam file with IGV to be sure
#fr-firststrand : R2 is in the same direction as the mRNA, and R1 is reverse
#fr-secondstrand : R1 is in the same direction as the mRNA, and R2 is reverse
#lib_type="fr-firststrand"

#bam list ; each of them will be subsampled (= a fraction a reads will be taken), and the results will be merged to have a metabam
#if there's a file in metabam_file variable, no need to fill this list, it will not be used
bam_list=()

#fraction of reads to take from each bam
#fract=0.2

#how many subprocesses should be run in parallel (subsampling of bam files) ?
#multiply this mumber by samtools_threads to have an idea of the max threads that will be used by this script
process_number_limit=1

#number of threads to give to samtools ; multiply this mumber by process_number_limit to have an idea of the max threads that will be used by this script
#samtools_threads=3

##############################################################

start_analysis=$(date)


############ -2)  process input variables #####################
##############################################################

#process input directory
output_dir="${output_dir}/"
output_dir=$(echo $output_dir |sed 's/\/\//\//g')
if [[ ! -d $output_dir ]]; then mkdir $output_dir; fi

#directory for subscripts
sample_scripts_dir="${output_dir}/sample_scripts/"
if [[ -d $sample_scripts_dir ]]; then rm -rf $sample_scripts_dir; fi
mkdir $sample_scripts_dir

if [[ ! -f $getFixedKmerStrand  ]];then

	getFixedKmerStrand=$(dirname "$0")/getFixedKmerStrand.awk

fi

if [[ ! -f $getFixedKmerStrand  ]];then

		echo -e "-> can't locate getFixedKmerStrand.awk script !\n"
		
		exit 1

fi


chmod 755 $getFixedKmerStrand

################################################################



######## -3)  create metabam if not supplied #######################
###################################################################

if [[ ! -f $metabam_file ]];then

		echo -e "we are going to create the metabam...\n"
		
		start_metabam=$(date)

		subbam_list=()
		
		if [[ ! -f ${output_dir}merged_sub_${fract}_bam ]];then
		

				for one_turn in $(seq 0 $((${#bam_list[*]}-1)));do
				
					if [[ ! -f ${output_dir}${one_turn}_subsample_${fract}_bam ]];then

						echo -e "$samtools view -@ $samtools_threads -bh -s $fract ${bam_list[$one_turn]} >${output_dir}${one_turn}_subsample_${fract}_bam" >${sample_scripts_dir}${one_turn}_subsampling_command.sh
				  
					fi
					
					subbam_list+=(${output_dir}${one_turn}_subsample_${fract}_bam)

				done
				
				find ${sample_scripts_dir} -name "*subsampling_*\.sh" | xargs -n 1 -P $process_number_limit bash || { echo "executing bash subscripts failure 1 !" 1>&2; exit; }
		  
				
				$samtools merge -@${samtools_threads} ${output_dir}merged_sub_${fract}_bam ${subbam_list[*]}
				
				
				end_metabam=$(date)
		
				echo -e "metabam created ! start : $start_metabam ; end : $end_metabam\n\n"
		  
		  

		else
		
				echo -e "metabam was already created (${output_dir}merged_sub_${fract}_bam) \n\n"
		
		
		fi
		
		#the metabam file is the merged bam from all subsampled bam files
		metabam_file=${output_dir}merged_sub_${fract}_bam

else

	echo -e "-> metabam already set to $metabam_file !\n"

fi

#####################################################################



############ -4)  split the metabam in plus strand bam and minus strand bam #######
##################################################################################


#build bam of plus strand
if [[ ! -f ${output_dir}merged_sub_${fract}_plus_bam  ]];then

	#fr-firststrand library
	##RNA
	#----------------------------------------->

	#R2 (flags -f128 -F16)					#R1 (flags -f64 -f16)
	#--->                                <----

	#see https://broadinstitute.github.io/picard/explain-flags.html for the flags

	#R2 forward strand, R1 reverse strand


    if [[ $lib_type == "fr-firststrand" ]];then

		#XS+
		echo -e "cat <($samtools view -h -@ $samtools_threads -f128 -F16 $metabam_file) <($samtools view -@ $samtools_threads -f64 -f16 $metabam_file) |$samtools sort -@ $samtools_threads -T ${output_dir}merged_sub_${fract}_bam_XS_plus_sorted_tmp - >${output_dir}merged_sub_${fract}_plus_bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >${sample_scripts_dir}subscript_2_1.sh
		
		
	#fr-secondstrand library
	##RNA
	#----------------------------------------->

	#R1 (flags -f64 -F16)					#R2 (flags -f128 -f16)
	#--->                                <----	
		
	#see https://broadinstitute.github.io/picard/explain-flags.html for the flags

	#R1 forward strand, R2 reverse strand
		
		
	
	elif [[ $lib_type == "fr-secondstrand" ]];then
	
		echo -e "cat <($samtools view -h -@ $samtools_threads -f128 -f16 $metabam_file) <($samtools view -@ $samtools_threads -f64 -F16 $metabam_file) |$samtools sort -@ $samtools_threads -T ${output_dir}merged_sub_${fract}_bam_XS_plus_sorted_tmp - >${output_dir}merged_sub_${fract}_plus_bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >${sample_scripts_dir}subscript_2_1.sh


	else
	
		echo -e "-> unknown lib_type $lib_type ! chose between fr-firststrand and fr-secondstrand\n"
		
		exit 1
	
	fi


else

	echo -e "${output_dir}merged_sub_${fract}_plus_bam  is already there\n\n"

fi

#build bam of minus strand
if [[ ! -f ${output_dir}merged_sub_${fract}_minus_bam  ]];then

	#fr-firststrand library
	##RNA
	##<-----------------------------------------

	##R1 (flags -f64 -F16)					#R2 (flags -f128 -f16)
	##--->                                <----

	#see https://broadinstitute.github.io/picard/explain-flags.html for the flags

	##R1 forward strand, R2 reverse strand

   if [[ $lib_type == "fr-firststrand" ]];then

		#XS-
		echo -e "cat <($samtools view -h -@ $samtools_threads -f128 -f16 $metabam_file) <($samtools view -@ $samtools_threads -f64 -F16 $metabam_file)  | $samtools sort -@ $samtools_threads -T ${output_dir}merged_sub_${fract}_bam_XS_minus_sorted_tmp - >${output_dir}merged_sub_${fract}_minus_bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >${sample_scripts_dir}subscript_2_2.sh
		
		
		
	#fr-secondstrand library
	##RNA
	##<-----------------------------------------

	##R2 (flags -f128 -F16)					#R1 (flags -f64 -f16)
	##--->                                <----

	#see https://broadinstitute.github.io/picard/explain-flags.html for the flags

	##R2 forward strand, R1 reverse strand		
		
	
	elif [[ $lib_type == "fr-secondstrand" ]];then
	
		#XS-
		echo -e "cat <($samtools view -h -@ $samtools_threads -f128 -F16 $metabam_file) <($samtools view -@ $samtools_threads -f64 -f16 $metabam_file)  | $samtools sort -@ $samtools_threads -T ${output_dir}merged_sub_${fract}_bam_XS_minus_sorted_tmp - >${output_dir}merged_sub_${fract}_minus_bam || { echo \"$samtools sort -@ $samtools_threads failure !\" 1>&2; exit; }\n" >${sample_scripts_dir}subscript_2_2.sh	
	
	
	else
	
		echo -e "-> unknown lib_type $lib_type ! chose between fr-firststrand and fr-secondstrand\n"
		
		exit 1	
	
	fi


else

	echo -e "${output_dir}merged_sub_${fract}_minus_bam  is already there\n\n"


fi


start_split=$(date)

find ${sample_scripts_dir} -name "*subscript_*\.sh" | xargs -n 1 -P 0 bash || { echo "executing bash subscripts failure 2 !" 1>&2; exit; }

end_split=$(date)

echo -e "metabam split strand + & - done ! start : $start_split ; end : $end_split\n\n"


#store just the mapped read sequence for plus strand to speed up the grep in downstream analysis
if [[ ! -f ${output_dir}merged_plus_seq.txt  ]];then

	$samtools view ${output_dir}merged_sub_${fract}_plus_bam|cut -f10 >${output_dir}merged_plus_seq.txt

fi


#store just the mapped read sequence for minus strand to speed up the grep in downstream analysis
if [[ ! -f ${output_dir}merged_minus_seq.txt  ]];then

	$samtools view ${output_dir}merged_sub_${fract}_minus_bam|cut -f10 >${output_dir}merged_minus_seq.txt


fi



#######################################################################################################################


############ -5)  use direct and indirect (= after reverse complement) grep on bam files (strand + & -)  to have counts ##########
#################################################################################################################################

#kmer seq (direct)
grep -v -i "contig" $kmer_file| cut -f${col_kmer} |LANG=en_EN sort -T ${output_dir} -k1,1 >${output_dir}query_kmers_direct.txt

#rev comp of kmer seq (indirect)
cat ${output_dir}query_kmers_direct.txt|rev|tr 'ATGC' 'TACG' >${output_dir}query_kmers_indirect.txt

#corresponding table rev comp of kmer seq and kmer seq
paste ${output_dir}query_kmers_indirect.txt ${output_dir}query_kmers_direct.txt |LANG=en_EN sort -T ${output_dir} -k1,1 >${output_dir}matching_kmer_seq_direct_and_indirect.tsv

#if the file with direct and indirect counts for both strand is empty, fill it
#put this if statement in comment if you want to rebuild the count files 
#if [[ ! -f ${output_dir}direct_and_indirect_kmer_counts.tsv ]] || [[ $(wc -l ${output_dir}direct_and_indirect_kmer_counts.tsv|awk '{print $1}') -lt $(wc -l ${output_dir}query_kmers_indirect.txt|awk '{print $1}')  ]];then

		grep_command_start=$(date)


		#grep candidate kmers (the commands will be run in parallel)
		#look for kmers on plus strand (direct)
		LC_ALL=C grep -o -F -f ${output_dir}query_kmers_direct.txt ${output_dir}merged_plus_seq.txt|LANG=en_EN sort -T ${output_dir}|uniq -c|awk 'OFS="\t"{print $2,$1}' >${output_dir}plus_strand_kmers_direct.tsv &


		#look for kmers on minus strand (direct)
		LC_ALL=C grep -o -F -f ${output_dir}query_kmers_direct.txt ${output_dir}merged_minus_seq.txt|LANG=en_EN sort -T ${output_dir}|uniq -c|awk 'OFS="\t"{print $2,$1}' >${output_dir}minus_strand_kmers_direct.tsv &


		#look for rev comp kmers on plus strand (indirect)
		LC_ALL=C grep -o -F -f ${output_dir}query_kmers_indirect.txt ${output_dir}merged_plus_seq.txt|LANG=en_EN sort -T ${output_dir}|uniq -c|awk 'OFS="\t"{print $2,$1}' >${output_dir}plus_strand_kmers_indirect.tsv &


		#look for rev comp kmers on minus strand (indirect)
		LC_ALL=C grep -o -F -f ${output_dir}query_kmers_indirect.txt ${output_dir}merged_minus_seq.txt|LANG=en_EN sort -T ${output_dir}|uniq -c|awk 'OFS="\t"{print $2,$1}' >${output_dir}minus_strand_kmers_indirect.tsv &

		#wait for the commands above
		wait

		#set 0 for missing plus kmers (direct)
		LANG=en_EN join -t $'\t' -a1 -11 -21 ${output_dir}query_kmers_direct.txt ${output_dir}plus_strand_kmers_direct.tsv|awk 'OFS="\t"{if($2==""){$2=0};print}' |LANG=en_EN sort -T ${output_dir} -k1,1 >${output_dir}plus_strand_kmers_direct.tmp && mv ${output_dir}plus_strand_kmers_direct.tmp ${output_dir}plus_strand_kmers_direct.tsv

		#set 0 for missing minus kmers (direct)
		LANG=en_EN join -t $'\t' -a1 -11 -21 ${output_dir}query_kmers_direct.txt ${output_dir}minus_strand_kmers_direct.tsv|awk 'OFS="\t"{if($2==""){$2=0};print}' |LANG=en_EN sort -T ${output_dir} -k1,1 >${output_dir}minus_strand_kmers_direct.tmp && mv ${output_dir}minus_strand_kmers_direct.tmp ${output_dir}minus_strand_kmers_direct.tsv

		#set 0 for missing plus kmers (indirect)
		LANG=en_EN join -t $'\t' -a1 -11 -21 ${output_dir}matching_kmer_seq_direct_and_indirect.tsv ${output_dir}plus_strand_kmers_indirect.tsv|awk 'OFS="\t"{if($3==""){$3=0};print $2,$3}' |LANG=en_EN sort -T ${output_dir} -k1,1 >${output_dir}plus_strand_kmers_indirect.tmp && mv ${output_dir}plus_strand_kmers_indirect.tmp ${output_dir}plus_strand_kmers_indirect.tsv

		#set 0 for missing minus kmers (indirect)
		LANG=en_EN join -t $'\t' -a1 -11 -21 ${output_dir}matching_kmer_seq_direct_and_indirect.tsv ${output_dir}minus_strand_kmers_indirect.tsv|awk 'OFS="\t"{if($3==""){$3=0};print $2,$3}' |LANG=en_EN sort -T ${output_dir} -k1,1 >${output_dir}minus_strand_kmers_indirect.tmp && mv ${output_dir}minus_strand_kmers_indirect.tmp ${output_dir}minus_strand_kmers_indirect.tsv

		#paste the count values in order to have in a table  <tag> <direct_plus_counts> <direct_minus_counts> <indirect_plus_counts> <indirect_minus_counts>
		paste ${output_dir}plus_strand_kmers_direct.tsv ${output_dir}minus_strand_kmers_direct.tsv ${output_dir}plus_strand_kmers_indirect.tsv ${output_dir}minus_strand_kmers_indirect.tsv|awk 'OFS="\t"{print $1,$2,$4,$6,$8}' >${output_dir}direct_and_indirect_kmer_counts.tsv

		grep_command_end=$(date)

		echo -e "start grep commands : $grep_command_start\t end grep commands : $grep_command_end"

#fi


############ -6)  determine the strandedness from direct and indirect counts, create fasta file accordingly ##########
###################################################################################################################


#corresponding table kmer <-> contig
paste <(grep -v -i "contig" $kmer_file| cut -f${col_kmer}) <(grep -v -i "contig" $kmer_file| cut -f${col_contig})|LANG=en_EN sort -T ${output_dir} -k1,1 >${output_dir}matching_kmer_contig.tsv


#add the contig info to the counts file
LANG=en_EN join -t $'\t' -a1 -11 -21 ${output_dir}direct_and_indirect_kmer_counts.tsv ${output_dir}matching_kmer_contig.tsv >${output_dir}direct_and_indirect_kmer_counts_annotated.tmp && mv ${output_dir}direct_and_indirect_kmer_counts_annotated.tmp ${output_dir}direct_and_indirect_kmer_counts_annotated.tsv

#determine the strandedness according to direct (strand + & -) and indirect (strand + & strand -) kmer counts
awk -f $getFixedKmerStrand ${output_dir}direct_and_indirect_kmer_counts_annotated.tsv >${output_dir}direct_and_indirect_kmer_counts_annotated.tmp && mv ${output_dir}direct_and_indirect_kmer_counts_annotated.tmp ${output_dir}direct_and_indirect_kmer_counts_annotated.tsv


#create fasta file
awk '{print ">"$8"\n"$9}' ${output_dir}direct_and_indirect_kmer_counts_annotated.tsv >${output_dir}direct_and_indirect_kmer_counts_annotated.fa

#add the header to the final table of counts
cat <(echo -e "tag\tplus_strand_direct_counts\tminus_strand_direct_counts\tplus_strand_indirect_counts\tminus_strand_indirect_counts\tselected_strand\toperation\tcontig_name\tselected_contig") ${output_dir}direct_and_indirect_kmer_counts_annotated.tsv >${output_dir}direct_and_indirect_kmer_counts_annotated.tmp && mv ${output_dir}direct_and_indirect_kmer_counts_annotated.tmp ${output_dir}direct_and_indirect_kmer_counts_annotated.tsv


####################################################################################################################

end_analysis=$(date)

echo -e "\n\ncheck files : \n${output_dir}direct_and_indirect_kmer_counts_annotated.tsv & ${output_dir}direct_and_indirect_kmer_counts_annotated.fa\n\nstart analysis : $start_analysis ; end analysis : $end_analysis"








