#!/usr/bin/env python
# coding: utf-8


import time
from Bio import SeqIO
from Bio import SearchIO
from collections import Counter

start_time=time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())
print("Primer check started at: %s\n"%start_time)

blast_qresults=SearchIO.parse('primerblast.xml','blast-xml',use_raw_hit_ids=True)

for blast_qresult in blast_qresults:
    primer_len=blast_qresult.seq_len
    seq_num=len(blast_qresult)
    mut_num=0
    low_freq_mut_num=0
    mut_list=[]
    
    for hit in blast_qresult:
        hsp=hit[0]
        if hsp.ident_num==primer_len:
            primer_seq=hsp.query.seq
            break
                
    for hit in blast_qresult:
        hsp=hit[0]
        if hsp.ident_num<primer_len:
                
            for reference in SeqIO.parse("gisaid_hcov-19.fasta", "fasta"):
                if hsp.hit_description in reference.id:
                    if hsp.hit_strand==1:
                        mut=reference.seq[hsp.hit_start-hsp.query_start:hsp.hit_start-hsp.query_start+primer_len].upper()
                    else:
                        mut=reference.seq[hsp.hit_start+hsp.query_end-primer_len:hsp.hit_start+hsp.query_end].reverse_complement().upper()
                    if mut.count("A")+mut.count("T")+mut.count("C")+mut.count("G")==primer_len:
                        mut_num+=1
                        mut_list.append(str(mut))
                    break
    
    mut_count=Counter(mut_list)
    
    print("%s\n%i mutants found in %i sequences\nP: %s"%(blast_qresult.id,mut_num,seq_num,primer_seq))
    for mut,freq in mut_count.most_common():
        if freq>1:
            print("M: %s F:%i"%(mut,freq))
        else:
            low_freq_mut_num+=1
    print("filtered %i low frequency mutations\n"%low_freq_mut_num)
    
end_time=time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())
print("Primer check finished at: %s"%end_time)

