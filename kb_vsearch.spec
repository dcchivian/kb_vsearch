/*
** A KBase module: kb_vsearch
**
** This module contains 4 methods from VSEARCH: basic query/db search, clustering, chimera detection, and dereplication.
** 
** Initially only basic query/db search will be implemented between read sets
*/

module kb_vsearch {
    typedef string workspace_id;

    /* 
    **    The workspace object refs are of form:
    **
    **    @id ws KBaseGenomes.ContigSet
    */
    typedef string ws_read_seq_set_ref; /* must be SingleEndLibrary */
    typedef string ws_feature_ref;
    typedef string ws_feature_set_ref;
    typedef string ws_genome_ref;
    typedef string ws_genome_set_ref;
    typedef string read_seq_set_name;
    typedef string feature_set_name;

    typedef string report_name;
    typedef string ws_report_ref;


    /* VSearch BasicSearch Input Params
    */
    typedef structure {
        workspace_id         workspace_id;
        ws_read_seq_set_ref  read_seq_set_one_ref;
        ws_read_seq_set_ref  read_seq_set_many_ref;
	ws_feature_ref       feature_one_ref;
	ws_feature_set_ref   feature_set_many_ref;
	ws_genome_ref        genome_many_ref;
	ws_genome_set_ref    genome_set_many_ref;

	read_seq_set_name    read_seq_set_output_name;
	feature_set_name     feature_set_output_name;

	int    maxaccepts;
	int    maxrejects;
	int    wordlength;
	int    minwordmatches;
	float  ident_thresh;
	int    ident_mode;
    } VSearch_BasicSearch_Params;


    /* VSearch BasicSearch Output
    */
    typedef structure {
	report_name    output_report_name;
	ws_report_ref  output_report_ref;

        ws_read_seq_set_ref      read_seq_set_output_ref;
        ws_feature_set_ref  feature_set_output_ref;

        int n_initial_seqs;
        int n_seqs_matched;
        int n_seqs_notmatched;
    } VSearch_BasicSearch_Output;
	
    /*
    **  Do basic search of one sequence against many sequences
    */
    funcdef VSearch_BasicSearch (VSearch_BasicSearch_Params params)  returns (VSearch_BasicSearch_Output) authentication required;
};
