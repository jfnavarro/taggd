

cdef object demultiplex_lines(str filename_reads,
                            str second_fastq_filename,
                            str filename_matched,
                            str filename_ambig,
                            str filename_unmatched,
                            str filename_res,
                            int ln_offset,
                            int ln_mod,
                            tuple umi_coordinates)

cdef list demultiplex_record(object rec)

cdef str trim_helpers(str seq)