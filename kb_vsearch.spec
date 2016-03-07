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
    **    objects = ws.get_objects([{'ref': params['workspace_id']+'/'+params['obj_name']}])
    **
    ** "ref" means the entire name combining the workspace id and the object name
    ** "id" is a numerical identifier of the workspace or object, and should just be used for workspace
    ** "name" is a string identifier of a workspace or object.  This is received from Narrative.
    */
    typedef string workspace_name;
    typedef string sequence;
    typedef string data_obj_name;
    typedef string data_obj_ref;


    /* VSearch BasicSearch Input Params
    */
    typedef structure {
        workspace_name workspace_name;
	sequence       input_one_string;
	data_obj_name  input_one_name;
	data_obj_name  input_many_name;
        data_obj_name  output_filtered_name;

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
	data_obj_name report_name;
	data_obj_ref  report_ref;
/*       data_obj_ref  output_filtered_ref;
*
*        int n_initial_seqs;
*        int n_seqs_matched;
*        int n_seqs_notmatched;
*/
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
