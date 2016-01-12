/*
** A KBase module: kb_vsearch
**
** This module contains 4 methods from VSEARCH: basic query/db search, clustering, chimera detection, and dereplication.
** 
** Initially only basic query/db search will be implemented between read sets
*/

module kb_vsearch {
    typedef string workspace_name;
    typedef string read_seq_set_id; /* must be SingleEndLibrary */
    typedef string feature_id;
    typedef string feature_set_id;
    typedef string genome_id;
    typedef string genome_set_id;

    typedef structure {
        workspace_name workspace;
        read_seq_set_id seq_set_one_id;
        read_seq_set_id seq_set_many_id;
	feature_id feature_one_id;
	feature_set_id feature_set_many_id;
	genome_id genome_many_id;
	genome_set_id genome_set_many_id;
    } VSearchBasicSearchParams;

    /* 
    **    The workspace ID for a ContigSet data object.
    **    @id ws KBaseGenomes.ContigSet
    */
    /*typedef string ws_contigset_id;*/

    typedef string ws_read_seq_set_id;
    /*typedef string ws_feature_id;*/   /* these do not exist independently? */
    typedef string ws_feature_set_id;
    typedef string ws_genome_id;
    typedef string ws_genome_set_id;


    typedef structure {
        ws_contigset_id new_contigset_ref;
        int n_initial_contigs;
        int n_seqs_below_threshold;
        int n_seqs_above_threshold;
    } VSearchBasicSearchResults
	
    /*
    **  Do basic search of one sequence against many sequences
    */
    funcdef vsearch_basic_search(VSearchBasicSearchParams params) returns (VSearchBasicSearchResults) authentication required;
};
