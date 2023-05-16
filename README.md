# get_coverage_from_bam 
# Need to add samtools software to environment variables

Contact: xul <xul@genepioneer.com>  <1275875706@qq.com>
Description:

        get reference seqence depth for every site and every base

Usage:
                perl /share/nas6/xul/program/chloroplast/RNAedit/script/get_coverage.pl \
                                -s test.bam \
                                -r reference.fa \
                                -p test \
                                -t DNA \
                                -o 01_run
  Options:

        -s --sam        <string>
                sam/bam file [must be sorted]

        -r --ref        <string>
                reference seqence file [fasta format]

        -p --prefix     <string>
                prefix

        -t --type       <string>
                input data type: DNA or RNA

        -q --qual       <int>
                Minimum mapping quality
                default: 40;

        -o --outdir     [string]
                default 01_run:
                        prefix_datatype_coverage.txt
                        prefix_datatype_mapped_consensus.txt
                        prefix_datatype_mapped_reads.txt

        -h      help
