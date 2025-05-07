#!/usr/bin/awk -f

#last update : 28/09/2023

#this script will be used by getStrandedKmers.sh ; it has to be in the same directory as this last

#function to reverse the sequence
function get_rev(my_string){					
					rev=""
					
			for (i=length(my_string);i>0;i--){
						
				rev=rev substr(my_string,i,1)
						
			}
			
			return rev
					
}

#function to complement the sequence
function get_comp(my_string){
	
	
    comp=gensub("A", "t","g",my_string); comp=gensub("T","a","g",comp); comp=gensub("G", "c","g",comp); comp=gensub("C", "g","g",comp); comp=toupper(comp);
    
    return comp
	
	
}


{

	    OFS="\t";
	    
	    counts_threshold=1;
	    
	    kmer=$1
	    
	    plus_strand_direct_counts=$2;
	    
	    minus_strand_direct_counts=$3;
	    
	    plus_strand_indirect_counts=$4;
	    
	    minus_strand_indirect_counts=$5;
	    
	    contig=$6;
	    
	    selected_contig=contig;
	    

        #strand from direct grep (no rev comp) ; initialization
		direct_strand=".";
		
		#strand from grep after rev comp ; initialization
		indirect_strand=".";

       #retained strand ; initialization 
		selected_strand=".";

        #nb reads associated with the strongest direct strand ; initialization
		strongest_direct_counts=0;

		#nb reads associated with the strongest rev comp strand ; initialization
		strongest_indirect_counts=0;

        #type of operation (direct from grep, or rev comp first) ; initialization
		operation=".";
				
		##############  compute variables from direct grep (without rev comp) #########

        #if plus strand reads from direct grep > 10*minus strand reads from direct grep
		if(plus_strand_direct_counts>(10*minus_strand_direct_counts) && plus_strand_direct_counts>=counts_threshold){
			
			#strongest direct counts has value of the plus strand
			strongest_direct_counts=plus_strand_direct_counts;
		
		    #direct strand is "+"
			direct_strand="+";
		
		#if minus strand reads from direct grep > 10*plus strand reads from direct grep
		}else if(minus_strand_direct_counts>(10*plus_strand_direct_counts) && minus_strand_direct_counts>=counts_threshold){
			
			#strongest direct counts has value of the minus strand
			strongest_direct_counts=minus_strand_direct_counts;
		
		    #direct strand is "-"
			direct_strand="-";
		
		}
		
		################################################################################
		
		
		##############  compute variables from grep after rev comp of the kmer #########
		
		#if plus strand reads from indirect grep > 10*minus strand reads from indirect grep
		if(plus_strand_indirect_counts>(10*minus_strand_indirect_counts) && plus_strand_indirect_counts>=counts_threshold){
			
			#strongest rev comp counts has value of the plus strand
			strongest_indirect_counts=plus_strand_indirect_counts;
			
		    indirect_strand="+";
		    
		#if minus strand reads from indirect grep > 10*plus strand reads from indirect grep
		}else if(minus_strand_indirect_counts>(10*plus_strand_indirect_counts) && minus_strand_indirect_counts>=counts_threshold){
			
			
			#strongest rev comp counts has value of the minus strand
			strongest_indirect_counts=minus_strand_indirect_counts;
			
			indirect_strand="-";
			
		}
		
		##################################################################################
		
		
		#decide between indirect and direct strand which one to chose
		#if the direct strand  (direct grep) or the strand from rev comp (indirect grep) aren't null, detremine which one gives the real strand
		if(direct_strand!="." || indirect_strand != "."){
			
			    #if value from direct grep > 10*value from indirect grep
				if(strongest_direct_counts>(10*strongest_indirect_counts)  && strongest_direct_counts>=counts_threshold){
				
				    #the retained strand is the one from the direct grep
					selected_strand=direct_strand;
				 
				    #gives the flag direct
					operation="direct";
				
				#if value from direct grep > 10*value from indirect grep
				}else if(strongest_indirect_counts>(10*strongest_direct_counts) && strongest_indirect_counts>=counts_threshold){
					
					#the retained strand is the one from the indirect grep
					selected_strand=indirect_strand;
					
					#gives the flag rev_comp
					operation="rev_comp";
					
					
					
				}
		 
		 }
		 
		 #if plus strand
		 if(selected_strand=="+"){
			 
			 contig_name=kmer"_pos";
			 
			 #if we had to rev comp the kmer, rev comp the contig
			 if(operation=="rev_comp"){
				 			 
				selected_contig=get_comp(get_rev(contig));
			 
		     }
		
		
		#if minus strand	 
		}else if(selected_strand=="-"){
			
			contig_name=kmer"_neg";
			
			#if we didn't have to rev comp the kmer, rev comp the contig
			if(operation!="rev_comp"){
						 
				selected_contig=get_comp(get_rev(contig));
			 
			 
		   }
		
		#if none of them, let unstranded	 
		}else{
			
			
			contig_name=kmer"_uns";
			
			
		}
		 

		 
		#print the kmer and infos
		print kmer,plus_strand_direct_counts,minus_strand_direct_counts,plus_strand_indirect_counts,minus_strand_indirect_counts,selected_strand,operation,contig_name,selected_contig;
	
}
