/*
** A KBase module: kb_vsearch
**
** This module contains 4 methods from VSEARCH: basic query/db search, clustering, chimera detection, and dereplication.
** 
** Initially only basic query/db search will be implemented between read sets
*/

module kb_vsearch {

    /* 
    ** The workspace object refs are of form:
    **
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_id']}])
    **
    ** "ref" means the entire name combining the workspace id and the object id
    ** "id" means just the portion of the ref corresponding to workspace or object
    ** "name" should not be used except as a field within an object
    */
    typedef string workspace_id;


    /* we will be overloading object types as follows:
    **
    **    input_one_id: SingleEndLibrary, FeatureSet
    **    input_many_id: SingleEndLibrary, FeatureSet, Genome, GenomeSet
    **    output_id: SingleEndLibrary (if input_many is SELib), FeatureSet
    */  
    typedef string one_id;
    typedef string many_id;
    typedef string output_id;
    typedef string output_ref;

    typedef string report_id;
    typedef string report_ref;


    /* VSearch BasicSearch Input Params
    */
    typedef structure {
        workspace_id  workspace_id;
	one_id        input_one_id;
	many_id       input_many_id;
        output_id     output_filtered_id;

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
	report_id   output_report_id;
	report_ref  output_report_ref;

        output_ref  output_filtered_ref;

        int n_initial_seqs;
        int n_seqs_matched;
        int n_seqs_notmatched;
    } VSearch_BasicSearch_Output;
	

    /*  Method for BasicSearch of one sequence against many sequences 
    **
    **    overloading as follows:
    **        input_one_id: SingleEndLibrary, FeatureSet
    **        input_many_id: SingleEndLibrary, FeatureSet, Genome, GenomeSet
    **        output_id: SingleEndLibrary (if input_many is SELib), FeatureSet
    */
    funcdef VSearch_BasicSearch (VSearch_BasicSearch_Params params)  returns (VSearch_BasicSearch_Output) authentication required;
};
